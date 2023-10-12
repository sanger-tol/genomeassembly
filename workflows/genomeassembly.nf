/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowGenomeassembly.initialise(params, log)

// TODO nf-core: Add all file path parameters for the pipeline to the list below
// Check input path parameters to see if they exist
def checkPathParamList = [ params.input ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }

if (params.bed_chunks_polishing) { bed_chunks_polishing = params.bed_chunks_polishing } else { bed_chunks_polishing = 100; }

if (params.cool_bin) { cool_bin = params.cool_bin } else { cool_bin = 1000; }

if (params.polishing_on) { polishing_on = params.polishing_on } else { polishing_on = false; }
if (params.hifiasm_hic_on) { hifiasm_hic_on = params.hifiasm_hic_on } else { hifiasm_hic_on = false; }
if ('organelles_on' in params.keySet() && !params.organelles_on) {  organelles_on = false } else { organelles_on = true; }
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { PREPARE_INPUT   } from '../subworkflows/local/prepare_input'
include { RAW_ASSEMBLY    } from '../subworkflows/local/raw_assembly' 
include { ORGANELLES            } from '../subworkflows/local/organelles' 
include { GENOMESCOPE_MODEL } from '../subworkflows/local/genomescope_model'
include { PURGE_DUPS as PURGE_DUPS_PRI      } from '../subworkflows/local/purge_dups'
include { PURGE_DUPS as PURGE_DUPS_ALT      } from '../subworkflows/local/purge_dups'
include { POLISHING       } from '../subworkflows/local/polishing'
include { SCAFFOLDING     } from '../subworkflows/local/scaffolding'
include { KEEP_SEQNAMES as KEEP_SEQNAMES_PRIMARY } from '../modules/local/keep_seqnames'
include { KEEP_SEQNAMES as KEEP_SEQNAMES_HAPLOTIGS } from '../modules/local/keep_seqnames'
include { HIC_MAPPING     } from '../subworkflows/local/hic_mapping'
include { GENOME_STATISTICS as GENOME_STATISTICS_RAW  } from '../subworkflows/local/assembly_stats'
include { GENOME_STATISTICS as GENOME_STATISTICS_RAW_HIC  } from '../subworkflows/local/assembly_stats'
include { GENOME_STATISTICS as GENOME_STATISTICS_PURGED  } from '../subworkflows/local/assembly_stats'
include { GENOME_STATISTICS as GENOME_STATISTICS_POLISHED  } from '../subworkflows/local/assembly_stats'
include { GENOME_STATISTICS as GENOME_STATISTICS_SCAFFOLDS } from '../subworkflows/local/assembly_stats'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { CAT_CAT as CAT_CAT_HAPLOTIGS } from "../modules/nf-core/cat/cat/main"
include { CAT_CAT as CAT_CAT_PURGEDUPS } from "../modules/nf-core/cat/cat/main"
include { SAMTOOLS_FAIDX as SAMTOOLS_FAIDX_PURGEDUPS   }  from '../modules/nf-core/samtools/faidx/main'

include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/custom/dumpsoftwareversions/main'
include { SEQTK_SUBSEQ as SEQTK_SUBSEQ_PRIMARY   } from '../modules/nf-core/seqtk/subseq/main'
include { SEQTK_SUBSEQ as SEQTK_SUBSEQ_HAPLOTIGS } from '../modules/nf-core/seqtk/subseq/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary

workflow GENOMEASSEMBLY {

    ch_versions = Channel.empty()

    //
    // SUBWORKFLOW: Read in yaml, validate and prepare for further steps
    //   
    PREPARE_INPUT(ch_input)
    ch_versions = ch_versions.mix(PREPARE_INPUT.out.versions)
        
    PREPARE_INPUT.out.hifi.set{ hifi_reads_ch }
 

    PREPARE_INPUT.out.hic.map{ meta, reads, motif -> reads }.set{ hic_reads_ch }

    GENOMESCOPE_MODEL( hifi_reads_ch )   

    RAW_ASSEMBLY( hifi_reads_ch , hic_reads_ch, hifiasm_hic_on )
    RAW_ASSEMBLY.out.primary_contigs.set{ primary_contigs_ch }
    RAW_ASSEMBLY.out.alternate_contigs.set{ haplotigs_ch }
    GENOME_STATISTICS_RAW( primary_contigs_ch.join(haplotigs_ch), 
                       PREPARE_INPUT.out.busco,
                       GENOMESCOPE_MODEL.out.hist,
                       GENOMESCOPE_MODEL.out.ktab
    )

    if ( hifiasm_hic_on ) {
        GENOME_STATISTICS_RAW_HIC( RAW_ASSEMBLY.out.primary_hic_contigs
                                    .join(RAW_ASSEMBLY.out.alternate_hic_contigs), 
                           PREPARE_INPUT.out.busco,
                           GENOMESCOPE_MODEL.out.hist,
                           GENOMESCOPE_MODEL.out.ktab
        )
    }
    hifi_reads_ch.join(primary_contigs_ch)
            .join(GENOMESCOPE_MODEL.out.model)
            .set{ purge_dups_input }
    PURGE_DUPS_PRI( purge_dups_input, 'primary' )
    PURGE_DUPS_PRI.out.pri.map{ meta, fasta -> [[id:meta.id], fasta] }
                          .set{ primary_contigs_ch }
    
    haplotigs_ch.combine( PURGE_DUPS_PRI.out.alt )
                    .map{ meta_h, h, meta_h_purged, h_purged -> [meta_h, [h, h_purged]]}
                    .set{ haplotigs_to_merge }
    
    CAT_CAT_HAPLOTIGS{ haplotigs_to_merge } 
    hifi_reads_ch.join(CAT_CAT_HAPLOTIGS.out.file_out)
            .join(GENOMESCOPE_MODEL.out.model)
            .set{ purge_dups_haploitgs_input }

    PURGE_DUPS_ALT( purge_dups_haploitgs_input, 'haplotigs' )

    PURGE_DUPS_ALT.out.pri.map{ meta, fasta -> [[id:meta.id], fasta] }
                        .set{ haplotigs_ch }
    GENOME_STATISTICS_PURGED( primary_contigs_ch.join(haplotigs_ch), 
                       PREPARE_INPUT.out.busco,
                       GENOMESCOPE_MODEL.out.hist,
                       GENOMESCOPE_MODEL.out.ktab
    )
    PURGE_DUPS_PRI.out.pri.combine(PURGE_DUPS_ALT.out.pri)
                        .map{ meta_pri, purged_pri, meta_alt, purged_alt -> [[id: meta_pri.id], [purged_pri, purged_alt]]}
                        .set{ purged_pri_alt_ch }
    CAT_CAT_PURGEDUPS( purged_pri_alt_ch )
    if ( organelles_on ) {
        if ( !polishing_on ) {
            ORGANELLES(hifi_reads_ch, CAT_CAT_PURGEDUPS.out.file_out, PREPARE_INPUT.out.mito)
        }
    }

    if ( polishing_on ) {
        SAMTOOLS_FAIDX_PURGEDUPS( CAT_CAT_PURGEDUPS.out.file_out, [[],[]] )
        CAT_CAT_PURGEDUPS.out.file_out.join( SAMTOOLS_FAIDX_PURGEDUPS.out.fai )
                                  .set{ reference_ch }

        PREPARE_INPUT.out.illumina_10X.map{ meta, reads, kmers -> [reads] }
                        .set{ illumina_10X_ch }
        
        POLISHING(reference_ch, illumina_10X_ch, bed_chunks_polishing)
        ch_versions = ch_versions.mix(POLISHING.out.versions)
        
        if ( organelles_on ) {
            ORGANELLES(hifi_reads_ch, POLISHING.out.fasta, PREPARE_INPUT.out.mito)
        }

        // Separate the primary and alternative contigs again after polishing
        // Separate primary contigs
        KEEP_SEQNAMES_PRIMARY(PURGE_DUPS_PRI.out.pri)
        ch_versions = ch_versions.mix(KEEP_SEQNAMES_PRIMARY.out.versions)
        SEQTK_SUBSEQ_PRIMARY(POLISHING.out.fasta, KEEP_SEQNAMES_PRIMARY.out.seqlist)
        ch_versions = ch_versions.mix(SEQTK_SUBSEQ_PRIMARY.out.versions)
        POLISHING.out.fasta.map{ meta, f -> [id: meta.id] }
                            .combine(SEQTK_SUBSEQ_PRIMARY.out.sequences)
                            .set{ primary_contigs_ch }
        
        // Separate alt contigs
        KEEP_SEQNAMES_HAPLOTIGS(PURGE_DUPS_ALT.out.pri)
        ch_versions = ch_versions.mix(KEEP_SEQNAMES_HAPLOTIGS.out.versions)
        SEQTK_SUBSEQ_HAPLOTIGS(POLISHING.out.fasta, KEEP_SEQNAMES_HAPLOTIGS.out.seqlist)
        ch_versions = ch_versions.mix(SEQTK_SUBSEQ_HAPLOTIGS.out.versions)
        POLISHING.out.fasta.map{ meta, f -> [id: meta.id] }
                            .combine(SEQTK_SUBSEQ_HAPLOTIGS.out.sequences)
                            .set{ haplotigs_contigs_ch }

        // Check genome stats for polished pri and alt
        primary_contigs_ch.join(haplotigs_contigs_ch)
                        .map{ meta, pri, alt -> [[id:meta.id], pri, alt]}
                        .set{ polished_asm_stats_input_ch }
        GENOME_STATISTICS_POLISHED( polished_asm_stats_input_ch, 
                       PREPARE_INPUT.out.busco,
                       GENOMESCOPE_MODEL.out.hist,
                       GENOMESCOPE_MODEL.out.ktab
        )
        ch_versions = ch_versions.mix(GENOME_STATISTICS_POLISHED.out.versions)
    }

    PREPARE_INPUT.out.hic.map{ meta, crams, motif -> [meta, crams] }
                         .set{ crams_ch }

    // Map HiC data to the primary assembly
    HIC_MAPPING ( primary_contigs_ch,crams_ch )
    ch_versions = ch_versions.mix(HIC_MAPPING.out.versions)

    SCAFFOLDING( HIC_MAPPING.out.bed, primary_contigs_ch, cool_bin )
    ch_versions = ch_versions.mix(SCAFFOLDING.out.versions)

    SCAFFOLDING.out.fasta.combine(haplotigs_ch)
                        .map{meta_s, fasta_s, meta_h, fasta_h -> [ meta_h, fasta_s, fasta_h ]}
                        .set{ stats_input_ch }

    GENOME_STATISTICS_SCAFFOLDS( stats_input_ch, 
                       PREPARE_INPUT.out.busco,
                       GENOMESCOPE_MODEL.out.hist,
                       GENOMESCOPE_MODEL.out.ktab
    )
    //
    // MODULE: Collate versions.yml file
    //
    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, null)
    }
    NfcoreTemplate.summary(workflow, params, log)
    if (params.hook_url) {
        NfcoreTemplate.IM_notification(workflow, params, summary_params, projectDir, log)
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
