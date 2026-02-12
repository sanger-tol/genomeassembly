include { GENOME_STATISTICS } from '../../../subworkflows/sanger-tol/genome_statistics'
include { HIFIASM_ASSEMBLY  } from '../../../subworkflows/local/hifiasm_assembly'
include { MITOHIFI_ASSEMBLY } from '../../../subworkflows/local/mitohifi_assembly'
include { PURGING           } from '../../../subworkflows/local/purging'
include { POLISHING         } from '../../../subworkflows/local/polishing'
include { SCAFFOLDING       } from '../../../subworkflows/local/scaffolding'

include { setupStage        } from '../../../functions/local/assembly_stages'

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
    // Logic: Separate out all the different unique specs, For the optional
    // stages (purging, polishing and scaffolding, organelles), we need to filter to
    // ensure that we have no null channels where these were unspecified
    //
    ch_output_specs = ch_specs.map { spec -> spec.subMap(["id", "hashes", "data"]) }.unique()
    ch_bin_assembly_specs = ch_specs.map { spec -> setupStage(spec, "bin_assembly") }.unique()
    ch_assembly_specs = ch_specs.map { spec -> setupStage(spec, "assembly") }.unique()
    ch_purging_specs = ch_specs.map { spec -> setupStage(spec, "purging") }.filter { spec -> spec }.unique()
    ch_polishing_specs = ch_specs.map { spec -> setupStage(spec, "polishing") }.filter { spec -> spec }.unique()
    ch_scaffolding_specs = ch_specs.map { spec -> setupStage(spec, "scaffolding") }.filter { spec -> spec }.unique()
    ch_organelle_specs = ch_specs.flatMap { spec -> [setupStage(spec, "mito"), setupStage(spec, "plastid")] }.filter { spec -> spec }.unique()

    //
    // Subworkflow: raw assembly of long reads using hifiasm
    //
    HIFIASM_ASSEMBLY(
        ch_bin_assembly_specs,
        ch_assembly_specs,
    )
    ch_assemblies = ch_assemblies.mix(HIFIASM_ASSEMBLY.out.hifiasm_assemblies)

    //
    // Subworkflow: purge assemblies with purge_dups pipeline
    //
    PURGING(
        ch_purging_specs,
        ch_assemblies,
        val_fastx_reads_per_chunk
    )
    ch_versions = ch_versions.mix(PURGING.out.versions)
    ch_assemblies = ch_assemblies.mix(PURGING.out.purged_assemblies)

    //
    // Subworkflow: polish assemblies with the polishing pipeline
    //
    POLISHING(
        ch_polishing_specs,
        ch_assemblies,
        val_polishing_container_provided,
        val_sequences_per_polishing_chunk
    )
    ch_versions = ch_versions.mix(POLISHING.out.versions)
    ch_assemblies = ch_assemblies.mix(POLISHING.out.polished_assemblies)

    //
    // Subworkflow: run hic-mapping and scaffolding
    //
    SCAFFOLDING(
        ch_scaffolding_specs,
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
            return [
                hash: spec.id,
                stage: spec.stage,
                data: spec.data,
                params: spec.params,
                stats: stats,
                busco: busco,
                merqury: merqury
            ]
        }

    //
    // Subworkflow: Find organellar genomes using mitohifi
    //
    MITOHIFI_ASSEMBLY(
        ch_organelle_specs,
        ch_assemblies,
    )

    //
    // Logic: Re-join stages to their input specifications for publishing
    //
    ch_hifiasm_out = ch_output_specs
        .combine(HIFIASM_ASSEMBLY.out.hifiasm_output)
        .filter { spec, hifiasm -> hifiasm.hash in spec.hashes.values() }
        .map { spec, hifiasm -> spec + hifiasm }

    ch_purging_out = ch_output_specs
        .combine(PURGING.out.purging_output)
        .filter { spec, purging -> purging.hash in spec.hashes.values() }
        .map { spec, purging -> spec + purging }

    ch_polishing_out = ch_output_specs
        .combine(POLISHING.out.polishing_output)
        .filter { spec, polishing -> polishing.hash in spec.hashes.values() }
        .map { spec, polishing -> spec + polishing }

    ch_scaffolding_out = ch_specs
        .combine(SCAFFOLDING.out.scaffolding_output)
        .filter { spec, scaffolding -> scaffolding.hash in spec.hashes.values() }
        .map { spec, scaffolding -> spec + scaffolding }

    ch_organelle_out = ch_specs
        .combine(MITOHIFI_ASSEMBLY.out.mitohifi_assemblies)
        .filter { spec, organelle -> organelle.hash in spec.hashes.values() }
        .map { spec, organelle -> spec + organelle }

    ch_statistics_out = ch_specs
        .combine(ch_statistics)
        .filter { spec, statistics -> statistics.hash in spec.hashes.values() }
        .map { spec, statistics -> spec + statistics }

    emit:
    hifiasm     = ch_hifiasm_out
    purging     = ch_purging_out
    polishing   = ch_polishing_out
    scaffolding = ch_scaffolding_out
    organelle   = ch_organelle_out
    statistics  = ch_statistics_out
    versions    = ch_versions

}
