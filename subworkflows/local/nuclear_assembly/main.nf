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
        .multiMap { spec, hap1, hap2, datasets ->
            def long_read_dataset = datasets.find { data ->
                data.id == spec.long_read_dataset && data.platform == spec.long_read_platform
            }
            def illumina_maternal_dataset = datasets.find { data ->
                data.id == spec.illumina_maternal_dataset && data.platform == "illumina"
            }
            def illumina_paternal_dataset = datasets.find { data ->
                data.id == spec.illumina_paternal_dataset && data.platform == "illumina"
            }
            def long_read_fastk = long_read_dataset.fastk
            def illumina_maternal_fastk = illumina_maternal_dataset?.fastk?.get(1) ?: []
            def illumina_paternal_fastk = illumina_paternal_dataset?.fastk?.get(1) ?: []

            assemblies: [meta, hap1, hap2]
            fastk: [meta, long_read_fastk.get(0), long_read_fastk.get(1), illumina_maternal_fastk, illumina_paternal_fastk]
            busco_lineage: [meta, meta.busco_lineage]
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
    // MITOHIFI(
    //     ch_organelle_specs,
    //     ch_assemblies,
    //     ch_data
    // )

    emit:
    versions = ch_versions

}
