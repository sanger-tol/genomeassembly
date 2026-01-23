#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    sanger-tol/genomeassembly
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/sanger-tol/genomeassembly
----------------------------------------------------------------------------------------
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { GENOMEASSEMBLY          } from './workflows/genomeassembly'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_genomeassembly_pipeline'
include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_nfcore_genomeassembly_pipeline'
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOWS FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Run main analysis pipeline depending on type of input
//
workflow SANGERTOL_GENOMEASSEMBLY {

    take:
    ch_specs
    ch_data
    val_kmer_size
    val_fastx_reads_per_chunk
    val_polishing_container_provided
    val_hic_aligner
    val_hic_mapping_cram_chunk_size
    val_scaffolding_cool_bin_size
    val_busco_lineage_directory

    main:

    //
    // WORKFLOW: Run pipeline
    //
    GENOMEASSEMBLY (
        ch_specs,
        ch_data,
        val_kmer_size,
        val_fastx_reads_per_chunk,
        val_polishing_container_provided,
        val_hic_aligner,
        val_hic_mapping_cram_chunk_size,
        val_scaffolding_cool_bin_size,
        val_busco_lineage_directory
    )
}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    main:
    //
    // SUBWORKFLOW: Run initialisation tasks
    //
    PIPELINE_INITIALISATION (
        params.version,
        params.validate_params,
        params.monochrome_logs,
        args,
        params.outdir,
        params.help,
        params.help_full,
        params.show_hidden,
        params.genomic_data,
        params.assembly_specs,
        params.polishing_longranger_container_path
    )

    //
    // WORKFLOW: Run main workflow
    //
    SANGERTOL_GENOMEASSEMBLY (
        PIPELINE_INITIALISATION.out.specs,
        PIPELINE_INITIALISATION.out.data,
        params.kmer_size,
        params.purging_reads_chunk_size,
        params.polishing_longranger_container_path,
        params.hic_aligner,
        params.hic_mapping_cram_chunk_size,
        params.scaffolding_cool_bin_size,
        params.busco_lineage_directory,
    )

    //
    // SUBWORKFLOW: Run completion tasks
    //
    PIPELINE_COMPLETION (
        params.email,
        params.email_on_fail,
        params.plaintext_email,
        params.outdir,
        params.monochrome_logs,
        params.hook_url,
    )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
