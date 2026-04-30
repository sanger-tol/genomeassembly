include { GENOME_STATISTICS } from '../../../subworkflows/sanger-tol/genome_statistics'
include { HIFIASM_ASSEMBLY  } from '../../../subworkflows/local/hifiasm_assembly'
include { MITOHIFI_ASSEMBLY } from '../../../subworkflows/local/mitohifi_assembly'
include { PURGING           } from '../../../subworkflows/local/purging'
include { POLISHING         } from '../../../subworkflows/local/polishing'
include { SCAFFOLDING       } from '../../../subworkflows/local/scaffolding'

include { setupStage        } from '../../../functions/local/assembly_stages'

include { GENERATE_SPECIFICATION_INDEX as INDEX_STAGE } from '../../../modules/local/generate_specification_index/main.nf'
include { GENERATE_SPECIFICATION_INDEX as INDEX_SPEC  } from '../../../modules/local/generate_specification_index/main.nf'

workflow NUCLEAR_ASSEMBLY {
    take:
    ch_specs
    val_fastx_reads_per_chunk
    val_polishing_container_provided
    val_sequences_per_polishing_chunk
    val_hic_aligner
    val_hic_mapping_cram_chunk_size
    val_scaffolding_cool_bin_size
    val_busco_lineage_directory

    main:
    ch_versions = channel.empty()
    ch_assemblies = channel.empty()

    //
    // Logic: Set up and deduplicate input stages
    //
    ch_stages = ch_specs.flatMap { spec ->
        spec.stages.collect { stageName, _stageData -> setupStage(spec, stageName) }
    }.unique()

    //
    // Subworkflow: raw assembly of long reads using hifiasm
    //
    HIFIASM_ASSEMBLY(
        ch_stages.filter { stage -> stage.stage == "base" },
        ch_stages.filter { stage -> stage.stage == "hifiasm_assembly" }
    )
    ch_assemblies = ch_assemblies.mix(HIFIASM_ASSEMBLY.out.hifiasm_assemblies)

    //
    // Subworkflow: purge assemblies with purge_dups pipeline
    //
    PURGING(
        ch_stages.filter { stage -> stage.stage == "purging" },
        ch_assemblies,
        val_fastx_reads_per_chunk
    )
    ch_assemblies = ch_assemblies.mix(PURGING.out.purged_assemblies)

    //
    // Subworkflow: polish assemblies with the polishing pipeline
    //
    POLISHING(
        ch_stages.filter { stage -> stage.stage == "polishing" },
        ch_assemblies,
        val_polishing_container_provided,
        val_sequences_per_polishing_chunk
    )
    ch_assemblies = ch_assemblies.mix(POLISHING.out.polished_assemblies)

    //
    // Subworkflow: run hic-mapping and scaffolding
    //
    SCAFFOLDING(
        ch_stages.filter { stage -> stage.stage == "scaffolding" },
        ch_assemblies,
        val_hic_aligner,
        val_hic_mapping_cram_chunk_size,
        val_scaffolding_cool_bin_size
    )
    ch_assemblies = ch_assemblies.mix(SCAFFOLDING.out.scaffolded_assemblies)

    //
    // Subworkflow: calculate genome statistcs
    //
    ch_genome_statistics_inputs = ch_assemblies
        .multiMap { spec, hap1, hap2 ->
            assemblies: [spec, hap1, hap2.size() > 0 ? hap2 : []]
            fastk: [spec, spec.data.long_read.fk_hist, spec.data.long_read.fk_ktab, spec.data.maternal.haptab, spec.data.paternal.haptab]
            busco_lineage: [spec, spec.params.busco_lineage]
        }

    GENOME_STATISTICS(
        ch_genome_statistics_inputs.assemblies,
        ch_genome_statistics_inputs.fastk,
        ch_genome_statistics_inputs.busco_lineage,
        val_busco_lineage_directory
    )

    ch_statistics = GENOME_STATISTICS.out.stats
        .join(GENOME_STATISTICS.out.busco, remainder: true)
        .join(GENOME_STATISTICS.out.merqury, remainder: true)
        .map { spec, stats, busco, merqury ->
            return spec.subMap(["id", "stage", "data", "params", "tools"]) + [
                hap: spec._hap,
                output: [
                    statistics: [
                        stats: stats,
                        busco: busco,
                        merqury: merqury
                    ]
                ]
            ]
        }

    //
    // Subworkflow: Find organellar genomes using mitohifi
    //
    MITOHIFI_ASSEMBLY(
        ch_stages.filter { stage -> stage.stage in ["mitohifi_mito", "mitohifi_plastid"] },
        ch_assemblies,
    )

    //
    // Module: Generate index files for each stage and the overall stage
    //
    ch_versions_to_index = channel.topic("versions")
        .map { task, tool, version -> [task.split(':').last(), tool, version] }
        .unique()
        .groupTuple(by: 0)
        .map { task, tools, versions ->
            [(task): [tools, versions].transpose().collectEntries()]
        }
        .reduce { a, b -> a + b }

    INDEX_STAGE(
        ch_stages.filter { stage -> stage.stage != "base" },
        ch_versions_to_index
    )

    INDEX_SPEC(
        ch_specs.map { spec -> spec.subMap(["name", "assembler", "data", "params", "tools"]) },
        ch_versions_to_index
    )

    //
    // Logic: Re-join stages to their input specifications for publishing
    //
    ch_output_specs = ch_specs.map { spec -> spec.subMap(["name", "hashes", "assembler"]) }.unique()

    ch_hifiasm_out = ch_output_specs
        .combine(HIFIASM_ASSEMBLY.out.hifiasm_output)
        .filter { spec, hifiasm -> hifiasm.id in spec.hashes.values() }
        .map { spec, hifiasm -> spec + hifiasm }

    ch_purging_out = ch_output_specs
        .combine(PURGING.out.purging_output)
        .filter { spec, purging -> purging.id in spec.hashes.values() }
        .map { spec, purging -> spec + purging }

    ch_polishing_out = ch_output_specs
        .combine(POLISHING.out.polishing_output)
        .filter { spec, polishing -> polishing.id in spec.hashes.values() }
        .map { spec, polishing -> spec + polishing }

    ch_scaffolding_out = ch_output_specs
        .combine(SCAFFOLDING.out.scaffolding_output)
        .filter { spec, scaffolding -> scaffolding.id in spec.hashes.values() }
        .map { spec, scaffolding -> spec + scaffolding }

    ch_organelle_out = ch_output_specs
        .combine(MITOHIFI_ASSEMBLY.out.mitohifi_assemblies)
        .filter { spec, organelle -> organelle.id in spec.hashes.values() }
        .map { spec, organelle -> spec + organelle }

    ch_statistics_out = ch_output_specs
        .combine(ch_statistics)
        .filter { spec, statistics -> statistics.id in spec.hashes.values() }
        .map { spec, statistics -> spec + statistics }

    ch_stage_indexes = ch_output_specs
        .combine(INDEX_STAGE.out.spec)
        .filter { spec, index -> index.id in spec.hashes.values() }
        .map { spec, index -> spec + index }

    emit:
    hifiasm       = ch_hifiasm_out
    purging       = ch_purging_out
    polishing     = ch_polishing_out
    scaffolding   = ch_scaffolding_out
    organelle     = ch_organelle_out
    statistics    = ch_statistics_out
    stage_indexes = ch_stage_indexes
    spec_indexes  = INDEX_SPEC.out.spec
    versions      = ch_versions
}
