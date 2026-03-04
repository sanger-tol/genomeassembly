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
include { specToAssemblyDir       } from './functions/local/publishing'
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
    stage_indexes    = GENOMEASSEMBLY.out.stage_indexes
    spec_indexes     = GENOMEASSEMBLY.out.spec_indexes
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
    stage_indexes    = SANGERTOL_GENOMEASSEMBLY.out.stage_indexes
    spec_indexes     = SANGERTOL_GENOMEASSEMBLY.out.spec_indexes
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
        path { spec ->
            spec.output.hifiasm.fasta >> "${spec.name}/raw/"
            spec.output.hifiasm.graphs >> "${spec.name}/raw/"
            spec.output.hifiasm.bed >> "${spec.name}/raw/"
            spec.output.hifiasm.log >> "${spec.name}/raw/"
        }
    }
    purging {
        path { spec ->
            spec.output.purging.fasta >> "${spec.name}/purging/"
            spec.output.purging.split_fa >> "${spec.name}/purging/split_aln/"
            spec.output.purging.split_paf >> "${spec.name}/purging/split_aln/"
            spec.output.purging.pbcstat_hist >> "${spec.name}/purging/coverage/"
            spec.output.purging.pbcstat_basecov >> "${spec.name}/purging/coverage/"
            spec.output.purging.calcuts_cutoffs >> "${spec.name}/purging/coverage/"
            spec.output.purging.calcuts_log >> "${spec.name}/purging/coverage/"
            spec.output.purging.histplot >> "${spec.name}/purging/coverage/"
            spec.output.purging.bed >> "${spec.name}/purging/purge_dups/"
            spec.output.purging.log >> "${spec.name}/purging/purge_dups/"
        }
    }
    polishing {
        path { spec ->
            spec.output.polishing.fasta >> "${spec.name}/polishing/"
            spec.output.polishing.longranger_bam >> "${spec.name}/polishing/aln/"
            spec.output.polishing.longranger_bai >> "${spec.name}/polishing/aln/"
            spec.output.polishing.longranger_csv >> "${spec.name}/polishing/aln/"
            spec.output.polishing.merged_vcf >> "${spec.name}/polishing/vcf/"
            spec.output.polishing.merged_vcf_tbi >> "${spec.name}/polishing/vcf/"
        }
    }
    scaffolding {
        path { spec ->
            spec.output.scaffolding.fasta >> "${spec.name}/scaffolding/"
            spec.output.scaffolding.bam >> "${spec.name}/scaffolding/hic_aln/"
            spec.output.scaffolding.bai >> "${spec.name}/scaffolding/hic_aln/"
            spec.output.scaffolding.stats >> "${spec.name}/scaffolding/hic_aln/"
            spec.output.scaffolding.flagstats >> "${spec.name}/scaffolding/hic_aln/"
            spec.output.scaffolding.idxstats >> "${spec.name}/scaffolding/hic_aln/"
            spec.output.scaffolding.yahs_agp >> "${spec.name}/scaffolding/yahs/"
            spec.output.scaffolding.yahs_bin >> "${spec.name}/scaffolding/yahs/"
            spec.output.scaffolding.yahs_initial >> "${spec.name}/scaffolding/yahs/"
            spec.output.scaffolding.yahs_intermeriate >> "${spec.name}/scaffolding/yahs/"
            spec.output.scaffolding.yahs_log >> "${spec.name}/scaffolding/yahs/"
            spec.output.scaffolding.pretext >> "${spec.name}/scaffolding/contact_maps/"
            spec.output.scaffolding.pretext_png >> "${spec.name}/scaffolding/contact_maps/"
            spec.output.scaffolding.cool >> "${spec.name}/scaffolding/contact_maps/"
            spec.output.scaffolding.hic >> "${spec.name}/scaffolding/contact_maps/"
        }
    }
    contigs_mitohifi {
        path { spec ->
            spec.output.mitohifi.mitohifi_reference_fa >> "${spec.name}/${spec.params.organelle}/ref/"
            spec.output.mitohifi.mitohifi_reference_gb >> "${spec.name}/${spec.params.organelle}/ref/"
            spec.output.mitohifi.mitohifi_files >> "${spec.name}/${spec.params.organelle}/"
        }
    }
    statistics {
        path { spec ->
            spec.output.statistics.stats >> [
                "${spec.name}",
                "${specToAssemblyDir(spec)}/"
            ].join("/")
            spec.output.statistics.merqury >> [
                "${spec.name}",
                "${specToAssemblyDir(spec)}",
                "merqury.${getPlatformShortName(spec.data.long_read.platform)}/"
            ].join("/")
            spec.output.statistics.busco >> [
                "${spec.name}",
                "${specToAssemblyDir(spec)}",
                "busco.${spec.params.busco_lineage}/"
            ].join("/")
        }
    }
    oatk {
        path { spec ->
            spec.output.oatk.mito_fasta >> "${spec.name}/"
            spec.output.oatk.pltd_fasta >> "${spec.name}/"
            spec.output.oatk.mito_bed >> "${spec.name}/"
            spec.output.oatk.pltd_bed >> "${spec.name}/"
            spec.output.oatk.mito_gfa >> "${spec.name}/"
            spec.output.oatk.pltd_gfa >> "${spec.name}/"
            spec.output.oatk.mito_annot >> "${spec.name}/"
            spec.output.oatk.pltd_annot >> "${spec.name}/"
            spec.output.oatk.initial_gfa >> "${spec.name}/"
            spec.output.oatk.final_gfa >> "${spec.name}/"
            spec.output.oatk.oatk_log >> "${spec.name}/"
            spec.output.oatk.bandage_plots >> "${spec.name}/"
        }
    }
    reads_mitohifi {
        path { spec ->
            spec.output.mitohifi.mitohifi_reference_fa >> "${spec.name}/ref/"
            spec.output.mitohifi.mitohifi_reference_gb >> "${spec.name}/ref/"
            spec.output.mitohifi.mitohifi_files >> "${spec.name}/"
        }
    }
    stage_indexes {
        path { spec ->
            spec.index >> "${spec.name}/${specToAssemblyDir(spec)}/index.json"
        }
    }
    spec_indexes {
        path { spec ->
            spec.index >> "${spec.name}/index.json"
        }
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
