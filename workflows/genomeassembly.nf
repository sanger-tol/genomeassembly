/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { PREPARE_INPUTS                            } from '../subworkflows/local/prepare_inputs'
include { NUCLEAR_ASSEMBLY                          } from '../subworkflows/local/nuclear_assembly'

// Functions
include { paramsSummaryMap                          } from 'plugin/nf-schema'
include { softwareVersionsToYAML                    } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText                    } from '../subworkflows/local/utils_nfcore_genomeassembly_pipeline'

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
    val_hic_slices_per_chunk
    val_busco_lineage
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
    ch_versions = ch_versions.mix(PREPARE_INPUTS.out.versions)

    //
    // Subworkflow: perform assembly of nuclear genome
    //
    NUCLEAR_ASSEMBLY(
        PREPARE_INPUTS.out.specs.filter { spec -> spec.assembler == "hifiasm" },
        PREPARE_INPUTS.out.data,
        PREPARE_INPUTS.out.merqury_trio_haptabs,
        val_busco_lineage,
        val_busco_lineage_directory,
        val_fastx_reads_per_chunk
    )
    ch_versions = ch_versions.mix(NUCLEAR_ASSEMBLY.out.versions)

    //
    // Collate and save software versions
    //
    def topic_versions = Channel.topic("versions")
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

    softwareVersionsToYAML(ch_versions.mix(topic_versions.versions_file))
        .mix(topic_versions_string)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name:  'genomeassembly_software_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }


    emit:
    versions       = ch_versions                 // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
