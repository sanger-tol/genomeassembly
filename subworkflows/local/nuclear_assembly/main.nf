include { GENOME_STATISTICS                         } from '../../../subworkflows/sanger-tol/genome_statistics'
include { HIFIASM_ASSEMBLY                          } from '../../../subworkflows/local/hifiasm_assembly'
include { PURGING                                   } from '../../../subworkflows/local/purging'
include { POLISHING                                 } from '../../../subworkflows/local/polishing'
include { SCAFFOLDING                               } from '../../../subworkflows/local/scaffolding'

workflow NUCLEAR_ASSEMBLY {
    take:
    ch_specs
    ch_data
    ch_merqury_trio_haptabs
    val_busco_lineage
    val_busco_lineage_directory
    val_fastx_reads_per_chunk

    main:
    ch_versions = channel.empty()
    ch_assemblies = channel.empty()

    //
    // Logic: Separate out all the different unique specs, For the optional
    // stages (purging, polishing and scaffolding), we need to filter to ensure
    // that we have no null channels where these were unspecified
    //
    ch_bin_assembly_specs = ch_specs.map { spec -> spec.bin_assembly }.unique()
    ch_assembly_specs = ch_specs.map { spec -> spec.assembly }.unique()
    ch_purging_specs = ch_specs.map { spec -> spec?.purging }.filter { spec -> spec }.unique()
    ch_polishing_specs = ch_specs.map { spec -> spec?.polishing }.filter { spec -> spec }.unique()
    ch_scaffolding_specs = ch_specs.map { spec -> spec?.scaffolding }.filter { spec -> spec }.unique()

    //
    // Subworkflow: raw assembly of long reads using hifiasm
    //
    HIFIASM_ASSEMBLY(
        ch_bin_assembly_specs,
        ch_assembly_specs,
        ch_data
    )
    ch_versions = ch_versions.mix(HIFIASM_ASSEMBLY.out.versions)
    ch_assemblies = ch_assemblies.mix(HIFIASM_ASSEMBLY.out.hifiasm_assemblies)

    //
    // Subworkflow: purge assemblies with purge_dups pipeline
    //
    PURGING(
        ch_purging_specs,
        ch_assemblies,
        ch_data,
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
        ch_data
    )
    ch_versions = ch_versions.mix(POLISHING.out.versions)
    ch_assemblies = ch_assemblies.mix(POLISHING.out.polished_assemblies)

    // //
    // // Logic: set up for scaffolding
    // //
    // ch_assemblies_for_scaffolding_split = ch_all_assemblies_after_polishing
    //     .branch { meta, _hap1, _hap2 ->
    //             def scaffold = false
    //             if(meta.assembly_stage == "polished") { scaffold = true }
    //             else if(meta.assembly_stage == "purged" && !meta.polish) { scaffold = true }
    //             else if(meta.assembly_stage == "raw" && !meta.purge && !meta.polish) { scaffold = true }

    //             scaffold: params.enable_scaffolding && scaffold
    //             no_scaffold: true
    //     }

    // //
    // // Subworkflow: run hic-mapping and scaffolding
    // //
    // SCAFFOLDING(
    //     ch_assemblies_for_scaffolding_split.scaffold,
    //     ch_hic_reads,
    //     params.hic_aligner,
    //     params.hic_mapping_cram_chunk_size,
    //     params.cool_bin
    // )
    // ch_versions = ch_versions.mix(SCAFFOLDING.out.versions)

    // //
    // // Logic: Mark scaffolded assemblies as scaffolded
    // //
    // ch_assemblies_scaffolded = SCAFFOLDING.out.assemblies
    //     .map { meta, asm1, asm2 ->
    //         def meta_new = meta + [assembly_stage: "scaffolded"]
    //         [meta_new, asm1, asm2]
    //     }

    // //
    // // Logic: collect all assemblies together
    // //
    // ch_all_assemblies = channel.empty()
    //     .mix(
    //         ch_assemblies_raw,
    //         ch_assemblies_purged,
    //         ch_assemblies_polished,
    //         ch_assemblies_scaffolded
    //     )

    // //
    // // Module: bgzip all output assemblies
    // //
    // ch_all_assemblies_to_bgzip = ch_all_assemblies
    //     .flatMap { meta, asm1, asm2 ->
    //         def meta_asm1 = meta + [_hap: "hap1"]
    //         def meta_asm2 = meta + [_hap: "hap2"]
    //         return [ [meta_asm1, asm1], [meta_asm2, asm2] ]
    //     }
    //     .filter { meta, asm -> asm }

    // BGZIP_ASSEMBLIES(ch_all_assemblies_to_bgzip)
    // ch_versions = ch_versions.mix(BGZIP_ASSEMBLIES.out.versions)

    // //
    // // Subworkflow: collect all assemblies and calculate assembly QC metrics
    // //
    // ch_hap1_for_statistics = BGZIP_ASSEMBLIES.out.output
    //     .flatMap { meta, asm ->
    //         meta._hap == "hap1" ? [[ meta - meta.subMap("_hap"), asm ]] : []
    //     }

    // ch_hap2_for_statisics = BGZIP_ASSEMBLIES.out.output
    //     .flatMap { meta, asm ->
    //         meta._hap == "hap2" ? [[ meta - meta.subMap("_hap"), asm ]] : []
    //     }

    // ch_assemblies_for_statistics = ch_hap1_for_statistics.join(ch_hap2_for_statisics, by: 0, remainder: true)

    // GENOME_STATISTICS(
    //     ch_assemblies_for_statistics,
    //     KMERS.out.fastk.collect(),
    //     KMERS.out.maternal_hapdb.collect(),
    //     KMERS.out.paternal_hapdb.collect(),
    //     val_busco_lineage,
    //     val_busco_lineage_directory
    // )
    // ch_versions = ch_versions.mix(GENOME_STATISTICS.out.versions)

    emit:
    versions = ch_versions

}
