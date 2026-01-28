include { GENOME_STATISTICS } from '../../../subworkflows/sanger-tol/genome_statistics'
include { HIFIASM_ASSEMBLY  } from '../../../subworkflows/local/hifiasm_assembly'
include { MITOHIFI_ASSEMBLY } from '../../../subworkflows/local/mitohifi_assembly'
include { PURGING           } from '../../../subworkflows/local/purging'
include { POLISHING         } from '../../../subworkflows/local/polishing'
include { SCAFFOLDING       } from '../../../subworkflows/local/scaffolding'

workflow NUCLEAR_ASSEMBLY {
    take:
    ch_specs
    val_fastx_reads_per_chunk
    val_polishing_container_provided
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
    ch_bin_assembly_specs = ch_specs.map { spec -> spec.stages.bin_assembly }.unique()
    ch_assembly_specs = ch_specs.map { spec -> spec.stages.assembly }.unique()
    ch_purging_specs = ch_specs.map { spec -> spec.stages?.purging }.filter { spec -> spec }.unique()
    ch_polishing_specs = ch_specs.map { spec -> spec.stages?.polishing }.filter { spec -> spec }.unique()
    ch_scaffolding_specs = ch_specs.map { spec -> spec.stages?.scaffolding }.filter { spec -> spec }.unique()
    ch_organelle_specs = ch_specs.flatMap { spec -> [spec.stages?.mito, spec.stages?.plastid] }.filter { spec -> spec }.unique()

    //
    // Subworkflow: raw assembly of long reads using hifiasm
    //
    HIFIASM_ASSEMBLY(
        ch_bin_assembly_specs,
        ch_assembly_specs,
    )
    ch_versions = ch_versions.mix(HIFIASM_ASSEMBLY.out.versions)
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
        val_polishing_container_provided
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
    ch_versions = ch_versions.mix(SCAFFOLDING.out.versions)
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
    ch_versions = ch_versions.mix(GENOME_STATISTICS.out.versions)

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
    ch_hifiasm_out = ch_specs.combine(HIFIASM_ASSEMBLY.out.hifiasm_output).combine(GENOME_STATISTICS.out.genome_statistics_output)
        .filter { meta, hifiasm, statistics ->
            hifiasm.id == statistics.id && hifiasm.id in meta.hashes.values()
        }
        .map { spec, hifiasm, statistics ->
            spec + [hifiasm: hifiasm, statistics: statistics]
        }.view()

    ch_purging_out = ch_specs.combine(PURGING.out.purging_output).combine(GENOME_STATISTICS.out.genome_statistics_output)
        .filter { meta, purging, statistics ->
            purging.id == statistics.id && purging.id in meta.hashes.values()
        }
        .map { spec, purging, statistics ->
            spec + [purging: purging, statistics: statistics]
        }

    ch_polishing_out = ch_specs.combine(POLISHING.out.polishing_output).combine(GENOME_STATISTICS.out.genome_statistics_output)
        .filter { meta, polishing, statistics ->
            purging.id == statistics.id && polishing.id in meta.hashes.values()
        }
        .map { spec, polishing, statistics ->
            spec + [polishing: polishing, statistics: statistics]
        }

    ch_scaffolding_out = ch_specs
        .combine(SCAFFOLDING.out.scaffolding_output.filter { res -> res._hap == "hap1" })
        .combine(SCAFFOLDING.out.scaffolding_output.filter { res -> res._hap == "hap2" })
        .combine(GENOME_STATISTICS.out.genome_statistics_output)
        .filter { meta, hap1, hap2, statistics ->
            (hap1.id == hap2.id) &&
            (hap1.id == statistics.id) &&
            hap1.id in meta.hashes.values()
        }
        .map { spec, hap1, hap2, statistics ->
            spec + [scaffolding: [hap1, hap2], statistics: statistics]
        }

    ch_organelle_out = ch_specs.combine(MITOHIFI_ASSEMBLY.out.mitohifi_assemblies)
        .filter { spec, mitohifi -> mitohifi.id in spec.hashes.values() }
        .map { spec, mitohifi ->
            spec + [mitohifi: mitohifi]
        }

    emit:
    hifiasm_output = ch_hifiasm_out
    versions = ch_versions

}
