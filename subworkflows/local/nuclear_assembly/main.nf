include { GENOME_STATISTICS                         } from '../../../subworkflows/sanger-tol/genome_statistics'
include { HIFIASM_ASSEMBLY                          } from '../../../subworkflows/local/hifiasm_assembly'
include { MITOHIFI                                  } from '../../../subworkflows/local/mitohifi'
include { PURGING                                   } from '../../../subworkflows/local/purging'
include { POLISHING                                 } from '../../../subworkflows/local/polishing'
include { SCAFFOLDING                               } from '../../../subworkflows/local/scaffolding'

workflow NUCLEAR_ASSEMBLY {
    take:
    ch_specs
    ch_data
    ch_merqury_trio_haptabs
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
    ch_bin_assembly_specs = ch_specs.map { spec -> spec.bin_assembly }.unique()
    ch_assembly_specs = ch_specs.map { spec -> spec.assembly }.unique()
    ch_purging_specs = ch_specs.map { spec -> spec?.purging }.filter { spec -> spec }.unique()
    ch_polishing_specs = ch_specs.map { spec -> spec?.polishing }.filter { spec -> spec }.unique()
    ch_scaffolding_specs = ch_specs.map { spec -> spec?.scaffolding }.filter { spec -> spec }.unique()
    ch_organelle_specs = ch_specs.flatMap { spec -> [spec?.mito, spec?.plastid] }.filter { spec -> spec }.unique()

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

    // Exclude polishing from the CLI output if not enabled to clean
    // up onscreen output
    if(val_polishing_container_provided) {
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
    }

    //
    // Subworkflow: run hic-mapping and scaffolding
    //
    SCAFFOLDING(
        ch_scaffolding_specs,
        ch_assemblies,
        ch_data,
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
        .combine(ch_data)
        .combine(ch_merqury_trio_haptabs.ifEmpty([[]]))
        .multiMap { spec, hap1, hap2, datasets, merqury_trio_haptabs ->
            def long_read_dataset = datasets.find { data ->
                data.id == spec.long_read_dataset && data.platform == spec.long_read_platform
            }
            def haptabs = merqury_trio_haptabs.find { data ->
                data.long_read_dataset == spec.long_read_dataset &&
                data.platform == spec.long_read_platform &&
                data.maternal_illumina_dataset == spec.illumina_maternal_dataset &&
                data.paternal_illumina_dataset == spec.illumina_paternal_dataset
            }

            def illumina_maternal_haptab = haptabs?.mat_haptab ?: []
            def illumina_paternal_haptab = haptabs?.pat_haptab ?: []

            def out_hap2 = hap2.size() > 0 ? hap2 : []

            assemblies: [spec, hap1, out_hap2]
            fastk: [spec, long_read_dataset.fk_hist, long_read_dataset.fk_ktab, illumina_maternal_haptab, illumina_paternal_haptab]
            busco_lineage: [spec, spec.busco_lineage]
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
    MITOHIFI(
        ch_organelle_specs,
        ch_assemblies,
        ch_data
    )

    emit:
    versions = ch_versions

}
