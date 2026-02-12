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

include { getPlatformShortName    } from './functions/local/publishing'
include { stageToAssemblyDir      } from './functions/local/publishing'
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
    val_sequences_per_polishing_chunk
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
        val_sequences_per_polishing_chunk,
        val_hic_aligner,
        val_hic_mapping_cram_chunk_size,
        val_scaffolding_cool_bin_size,
        val_busco_lineage_directory
    )

    emit:
    datasets         = GENOMEASSEMBLY.out.datasets
    hifiasm          = GENOMEASSEMBLY.out.hifiasm
    purging          = GENOMEASSEMBLY.out.purging
    polishing        = GENOMEASSEMBLY.out.polishing
    scaffolding      = GENOMEASSEMBLY.out.scaffolding
    contigs_mitohifi = GENOMEASSEMBLY.out.contigs_mitohifi
    oatk             = GENOMEASSEMBLY.out.oatk
    reads_mitohifi   = GENOMEASSEMBLY.out.reads_mitohifi
    statistics       = GENOMEASSEMBLY.out.statistics
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
        params.polishing_n_chunks_bed,
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

    publish:
    datasets         = SANGERTOL_GENOMEASSEMBLY.out.datasets
    hifiasm          = SANGERTOL_GENOMEASSEMBLY.out.hifiasm
    purging          = SANGERTOL_GENOMEASSEMBLY.out.purging
    polishing        = SANGERTOL_GENOMEASSEMBLY.out.polishing
    scaffolding      = SANGERTOL_GENOMEASSEMBLY.out.scaffolding
    statistics       = SANGERTOL_GENOMEASSEMBLY.out.statistics
    contigs_mitohifi = SANGERTOL_GENOMEASSEMBLY.out.contigs_mitohifi
    oatk             = SANGERTOL_GENOMEASSEMBLY.out.oatk
    reads_mitohifi   = SANGERTOL_GENOMEASSEMBLY.out.reads_mitohifi
}

output {
    datasets {
        path { dataset ->
            dataset.fk_hist >> "fastk/${dataset.id}/${dataset.platform}/k${params.kmer_size}/"
            dataset.fk_ktab >> "fastk/${dataset.id}/${dataset.platform}/k${params.kmer_size}/"
        }
        enabled params.save_fastk_databases
    }
    hifiasm {
        path { assembly ->
            assembly.fasta >> "${assembly.id}/raw/"
            assembly.graphs >> "${assembly.id}/raw/"
            assembly.bed >> "${assembly.id}/raw/"
            assembly.log >> "${assembly.id}/raw/"
        }
    }
    purging {
        path { assembly ->
            assembly.fasta >> "${assembly.id}/purging/"
            assembly.split_fa >> "${assembly.id}/purging/split_aln/"
            assembly.split_paf >> "${assembly.id}/purging/split_aln/"
            assembly.pbcstat_hist >> "${assembly.id}/purging/coverage/"
            assembly.pbcstat_basecov >> "${assembly.id}/purging/coverage/"
            assembly.calcuts_cutoffs >> "${assembly.id}/purging/coverage/"
            assembly.calcuts_log >> "${assembly.id}/purging/coverage/"
            assembly.histplot >> "${assembly.id}/purging/coverage/"
            assembly.bed >> "${assembly.id}/purging/purge_dups/"
            assembly.log >> "${assembly.id}/purging/purge_dups/"
        }
    }
    polishing {
        path { assembly ->
            assembly.fasta >> "${assembly.id}/polishing/"
            assembly.bed_chunks >> "${assembly.id}/polishing/chunks/"
            assembly.longranger_bam >> "${assembly.id}/polishing/aln/"
            assembly.longranger_bai >> "${assembly.id}/polishing/aln/"
            assembly.longranger_csv >> "${assembly.id}/polishing/aln/"
            assembly.merged_vcf >> "${assembly.id}/polishing/vcf/"
            assembly.merged_vcf_tbi >> "${assembly.id}/polishing/vcf/"
        }
    }
    scaffolding {
        path { assembly ->
            assembly.fasta >> "${assembly.id}/scaffolding/${assembly.hap}/yahs/"
            assembly.bam >> "${assembly.id}/scaffolding/${assembly.hap}/hic_aln/"
            assembly.bai >> "${assembly.id}/scaffolding/${assembly.hap}/hic_aln/"
            assembly.stats >> "${assembly.id}/scaffolding/${assembly.hap}/hic_aln/"
            assembly.flagstats >> "${assembly.id}/scaffolding/${assembly.hap}/hic_aln/"
            assembly.idxstats >> "${assembly.id}/scaffolding/${assembly.hap}/hic_aln/"
            assembly.yahs_agp >> "${assembly.id}/scaffolding/${assembly.hap}/yahs/"
            assembly.yahs_bin >> "${assembly.id}/scaffolding/${assembly.hap}/yahs/"
            assembly.yahs_initial >> "${assembly.id}/scaffolding/${assembly.hap}/yahs/"
            assembly.yahs_intermeriate >> "${assembly.id}/scaffolding/${assembly.hap}/yahs/"
            assembly.yahs_log >> "${assembly.id}/scaffolding/${assembly.hap}/yahs/"
            assembly.pretext >> "${assembly.id}/scaffolding/${assembly.hap}/yahs/"
            assembly.pretext_png >> "${assembly.id}/scaffolding/${assembly.hap}/yahs/"
            assembly.cool >> "${assembly.id}/scaffolding/${assembly.hap}/yahs/"
            assembly.hic >> "${assembly.id}/scaffolding/${assembly.hap}/yahs/"
        }
    }
    contigs_mitohifi {
        path { assembly ->
            assembly.mitohifi_reference_fa >> "${assembly.id}/${assembly.params.organelle}/ref/"
            assembly.mitohifi_reference_gb >> "${assembly.id}/${assembly.params.organelle}/ref/"
            assembly.mitohifi_files >> "${assembly.id}/${assembly.params.organelle}/"
        }
    }
    statistics {
        path { statistics ->
            statistics.stats >> [
                "${statistics.id}",
                "${stageToAssemblyDir(statistics.stage)}/"
            ].join("/")
            statistics.merqury >> [
                "${statistics.id}",
                "${stageToAssemblyDir(statistics.stage)}",
                "merqury.${getPlatformShortName(statistics.data.long_read.platform)}/"
            ].join("/")
            statistics.busco >> [
                "${statistics.id}",
                "${stageToAssemblyDir(statistics.stage)}",
                "busco.${statistics.params.busco_lineage}/"
            ].join("/")
        }
    }
    oatk {
        path { assembly ->
            assembly.mito_fasta >> "${assembly.id}/"
            assembly.pltd_fasta >> "${assembly.id}/"
            assembly.mito_bed >> "${assembly.id}/"
            assembly.pltd_bed >> "${assembly.id}/"
            assembly.mito_gfa >> "${assembly.id}/"
            assembly.pltd_gfa >> "${assembly.id}/"
            assembly.mito_annot >> "${assembly.id}/"
            assembly.pltd_annot >> "${assembly.id}/"
            assembly.initial_gfa >> "${assembly.id}/"
            assembly.final_gfa >> "${assembly.id}/"
            assembly.oatk_log >> "${assembly.id}/"
            assembly.bandage_plots >> "${assembly.id}/"
        }
    }
    reads_mitohifi {
        path { assembly ->
            assembly.mitohifi_reference_fa >> "${assembly.id}/ref/"
            assembly.mitohifi_reference_gb >> "${assembly.id}/ref/"
            assembly.mitohifi_files >> "${assembly.id}/"
        }
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
