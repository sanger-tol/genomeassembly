/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { PREPARE_INPUTS         } from '../subworkflows/local/prepare_inputs'
include { NUCLEAR_ASSEMBLY       } from '../subworkflows/local/nuclear_assembly'
include { ORGANELLE_ASSEMBLY     } from '../subworkflows/local/organelle_assembly'

// Functions
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_genomeassembly_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow GENOMEASSEMBLY {

    take:
    ch_specs
    ch_data
    val_kmer_size
    val_fastx_reads_per_chunk
    val_polishing_container_provided
    val_sequences_per_polishing_chunk
    val_hic_aligner
    val_hic_mapping_cram_chunk_size
    val_scaffolding_cool_bin_size
    val_busco_lineage_directory

    main:
    ch_versions = channel.empty()

    //
    // Subworkflow: build kmer databases for input data if required
    // and then pre-process the specifications for downstream consumption
    //
    PREPARE_INPUTS(
        ch_specs,
        ch_data,
        val_kmer_size
    )

    //
    // Subworkflow: perform assembly of nuclear genome
    //
    NUCLEAR_ASSEMBLY(
        PREPARE_INPUTS.out.specs.filter { spec -> spec.assembler in ["hifiasm"] },
        val_fastx_reads_per_chunk,
        val_polishing_container_provided,
        val_sequences_per_polishing_chunk,
        val_hic_aligner,
        val_hic_mapping_cram_chunk_size,
        val_scaffolding_cool_bin_size,
        val_busco_lineage_directory
    )
    ch_versions = ch_versions.mix(NUCLEAR_ASSEMBLY.out.versions)

    //
    // Subworkflow: assemble organellar genomes
    //
    ORGANELLE_ASSEMBLY(
        PREPARE_INPUTS.out.specs.filter { spec -> spec.assembler in ["oatk", "mitohifi"] }
    )

    //
    // Collate and save software versions
    //
    def topic_versions = channel.topic("versions")
        .distinct()
        .branch { entry ->
            versions_file: entry instanceof Path
            versions_tuple: true
        }

    def topic_versions_string = topic_versions.versions_tuple
        .map { process, tool, version ->
            [ process[process.lastIndexOf(':')+1..-1], "  ${tool}: ${version}" ]
        }
        .groupTuple(by:0)
        .map { process, tool_versions ->
            tool_versions.unique().sort()
            "${process}:\n${tool_versions.join('\n')}"
        }

    ch_collated_versions = softwareVersionsToYAML(ch_versions.mix(topic_versions.versions_file))
        .mix(topic_versions_string)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name:  'genomeassembly_software_'  + 'versions.yml',
            sort: true,
            newLine: true
        )

    emit:
    datasets         = PREPARE_INPUTS.out.datasets
    hifiasm          = NUCLEAR_ASSEMBLY.out.hifiasm
    purging          = NUCLEAR_ASSEMBLY.out.purging
    polishing        = NUCLEAR_ASSEMBLY.out.polishing
    scaffolding      = NUCLEAR_ASSEMBLY.out.scaffolding
    contigs_mitohifi = NUCLEAR_ASSEMBLY.out.organelle
    statistics       = NUCLEAR_ASSEMBLY.out.statistics
    oatk             = ORGANELLE_ASSEMBLY.out.oatk
    reads_mitohifi   = ORGANELLE_ASSEMBLY.out.reads_mitohifi
    stage_indexes    = NUCLEAR_ASSEMBLY.out.stage_indexes
    spec_indexes     = NUCLEAR_ASSEMBLY.out.spec_indexes.mix(ORGANELLE_ASSEMBLY.out.spec_indexes)
    versions         = ch_collated_versions // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
