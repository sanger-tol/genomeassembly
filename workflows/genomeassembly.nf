/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Check input path parameters to see if they exist
def checkPathParamList = [ params.input ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input YAML not specified!' }

if (params.bed_chunks_polishing) { bed_chunks_polishing = params.bed_chunks_polishing } else { bed_chunks_polishing = 100; }

if (params.cool_bin) { cool_bin = params.cool_bin } else { cool_bin = 1000; }

if (params.polishing_on) { polishing_on = params.polishing_on } else { polishing_on = false; }
if (params.hifiasm_hic_on) { hifiasm_hic_on = params.hifiasm_hic_on } else { hifiasm_hic_on = false; }
if (params.organelles_on) { organelles_on = params.organelles_on } else { organelles_on = false; }
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
include { PREPARE_INPUT                                    } from '../subworkflows/local/prepare_input'
include { RAW_ASSEMBLY                                     } from '../subworkflows/local/raw_assembly' 
include { ORGANELLES                                       } from '../subworkflows/local/organelles' 
include { GENOMESCOPE_MODEL                                } from '../subworkflows/local/genomescope_model'
include { PURGE_DUPS                                       } from '../subworkflows/local/purge_dups'
include { POLISHING                                        } from '../subworkflows/local/polishing'
include { SCAFFOLDING                                      } from '../subworkflows/local/scaffolding'
include { KEEP_SEQNAMES as KEEP_SEQNAMES_PRIMARY           } from '../modules/local/keep_seqnames'
include { KEEP_SEQNAMES as KEEP_SEQNAMES_HAPLOTIGS         } from '../modules/local/keep_seqnames'
include { HIC_MAPPING                                      } from '../subworkflows/local/hic_mapping'
include { GENOME_STATISTICS as GENOME_STATISTICS_RAW       } from '../subworkflows/local/genome_statistics'
include { GENOME_STATISTICS as GENOME_STATISTICS_RAW_HIC   } from '../subworkflows/local/genome_statistics'
include { GENOME_STATISTICS as GENOME_STATISTICS_PURGED    } from '../subworkflows/local/genome_statistics'
include { GENOME_STATISTICS as GENOME_STATISTICS_POLISHED  } from '../subworkflows/local/genome_statistics'
include { GENOME_STATISTICS as GENOME_STATISTICS_SCAFFOLDS } from '../subworkflows/local/genome_statistics'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { CAT_CAT as CAT_CAT_MITOHIFI_READS          } from "../modules/nf-core/cat/cat/main"
include { CAT_CAT as CAT_CAT_RAW                     } from "../modules/nf-core/cat/cat/main"
include { CAT_CAT as CAT_CAT_HAPLOTIGS               } from "../modules/nf-core/cat/cat/main"
include { CAT_CAT as CAT_CAT_PURGEDUPS               } from "../modules/nf-core/cat/cat/main"
include { SAMTOOLS_FAIDX as SAMTOOLS_FAIDX_PURGEDUPS } from '../modules/nf-core/samtools/faidx/main'
include { SEQTK_SUBSEQ as SEQTK_SUBSEQ_HAPLOTIGS     } from '../modules/nf-core/seqtk/subseq/main'
include { SEQTK_SUBSEQ as SEQTK_SUBSEQ_PRIMARY       } from '../modules/nf-core/seqtk/subseq/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS                } from '../modules/nf-core/custom/dumpsoftwareversions/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary

workflow GENOMEASSEMBLY {

    ch_versions = Channel.empty()

    //
    // SUBWORKFLOW: READ IN YAML, VALIDATE AND PREPARE FOR FURTHER STEPS
    //   
    PREPARE_INPUT(ch_input)
    ch_versions = ch_versions.mix(PREPARE_INPUT.out.versions)
    
    //
    // LOGIC: CREATE A VARIABLE SERVING AS AN ALIAS FOR HIFI READS CHANNEL
    //
    PREPARE_INPUT.out.hifi.set{ hifi_reads_ch }
    
    //
    // LOGIC: SEPARATE READS PATHS INTO A DIFFERENT CHANNEL
    //    
    PREPARE_INPUT.out.hic.map{ meta, reads, motif -> reads }.set{ hic_reads_ch }

    //
    // SUBWORKFLOW: GENERATE KMER DATABASE AND PROFILE MODEL
    //
    GENOMESCOPE_MODEL( hifi_reads_ch )   
    ch_versions = ch_versions.mix(GENOMESCOPE_MODEL.out.versions)

    //
    // SUBWORKFLOW: RUN A HIFIASM ASSEMBLY ON THE HIFI READS; ALSO CREATE
    //              A HIFIASM RUN IN HIC MODE IF THE FLAG IS SWITCHED ON
    //
    RAW_ASSEMBLY( hifi_reads_ch, hic_reads_ch, hifiasm_hic_on )
    ch_versions = ch_versions.mix(RAW_ASSEMBLY.out.versions)

    //
    // LOGIC: DEFINE THE PRIMARY CONTIGS CHANNEL
    //
    RAW_ASSEMBLY.out.primary_contigs.set{ primary_contigs_ch }
    
    //
    // LOGIC: DEFINE THE HAPLOTIGS CHANNELS
    //
    RAW_ASSEMBLY.out.alternate_contigs.set{ haplotigs_ch }

    //
    // SUBWORKFLOW: CALCULATE STATISTICS FOR THE RAW ASSEMBLY
    //
    GENOME_STATISTICS_RAW( primary_contigs_ch.join(haplotigs_ch), 
                       PREPARE_INPUT.out.busco,
                       GENOMESCOPE_MODEL.out.hist,
                       GENOMESCOPE_MODEL.out.ktab
    )
    ch_versions = ch_versions.mix(GENOME_STATISTICS_RAW.out.versions)

    if ( organelles_on ) {
        //
        // LOGIC: CREATE CHANNEL FOR PRIMARY AND ALT CONTIGS
        //
        primary_contigs_ch.join(haplotigs_ch)
                          .map{ meta, pri, alt -> [meta, [pri, alt]]}
                          .set{ raw_pri_alt_ch }
        //
        // MODULE: MERGE PAW CONTIGS AND HAPLOTIGS INTO ONE FILE
        //
        CAT_CAT_RAW( raw_pri_alt_ch )

        //
        // LOGIC: DEFINE MERGED ASSEMBLY
        //
        merged_pri_alt_raw = CAT_CAT_RAW.out.file_out

        //
        // MODULE: MERGE INPUT FASTA FILES WITH PACBIO READS
        //
        CAT_CAT_MITOHIFI_READS(hifi_reads_ch)
        ch_versions = ch_versions.mix(CAT_CAT_MITOHIFI_READS.out.versions)

        //
        // SUBWORKFLOW: INDETIFY MITO IN THE RAW READS AND ASSEMBLY CONTIGS
        // 
        ORGANELLES(CAT_CAT_MITOHIFI_READS.out.file_out, merged_pri_alt_raw, PREPARE_INPUT.out.mito)
    }

    //
    // LOGIC: CHECK IF THE HIFIASM HIC MODE WAS SWITCHED ON
    //
    if ( hifiasm_hic_on ) {
        //
        // SUBWORKFLOW: CALCULATE RAW ASSEMBLY STATISTICS FOR THE HIFIASN IN HIC MODE
        //
        GENOME_STATISTICS_RAW_HIC( RAW_ASSEMBLY.out.primary_hic_contigs
                                    .join(RAW_ASSEMBLY.out.alternate_hic_contigs), 
                           PREPARE_INPUT.out.busco,
                           GENOMESCOPE_MODEL.out.hist,
                           GENOMESCOPE_MODEL.out.ktab
        )
    }

    //
    // LOGIC: CREATE AN INPUT DATA STRUCTURE FOR PURGING
    //
    hifi_reads_ch.join(primary_contigs_ch)
            .join(GENOMESCOPE_MODEL.out.model)
            .set{ purge_dups_input }

    //
    // SUBWORKFLOW: RUN PURGE DUPS ON THE PRIMARY CONTIGS
    //
    PURGE_DUPS( purge_dups_input )
    ch_versions = ch_versions.mix(PURGE_DUPS.out.versions)

    //
    // LOGIC: UPDATE THE PRIMARY CONTIGS CHANNEL
    //
    PURGE_DUPS.out.pri.map{ meta, fasta -> [[id:meta.id], fasta] }
                          .set{ primary_contigs_ch }
    
    //
    // LOGIC: SET APART THE HAPLOTIGS AFTER PURGING AND THE HIFIASM HAPLOTIGS
    //
    haplotigs_ch.combine( PURGE_DUPS.out.alt )
                    .map{ meta_h, h, meta_h_purged, h_purged -> [meta_h, [h, h_purged]]}
                    .set{ haplotigs_to_merge }
    
    //
    // MODULE: COMBINE PURGED SEQUENCES WITH THE ORIGINAL HAPLOTIGS
    //
    CAT_CAT_HAPLOTIGS{ haplotigs_to_merge }
    haplotigs_ch = CAT_CAT_HAPLOTIGS.out.file_out

    //
    // SUBWORKFLOW: CALCULATE STATISTICS FOR THE PURGED ASSEMBLY
    //
    GENOME_STATISTICS_PURGED( primary_contigs_ch.join(haplotigs_ch), 
                       PREPARE_INPUT.out.busco,
                       GENOMESCOPE_MODEL.out.hist,
                       GENOMESCOPE_MODEL.out.ktab
    )
   
    //
    // LOGIC: CREATE A CHANNEL FOR THE PURGED CONTIGS AMD HAPLOTIGS 
    //
    PURGE_DUPS.out.pri.join(haplotigs_ch)
                        .map{ meta, purged_pri, purged_alt -> [meta, [purged_pri, purged_alt]]}
                        .set{ purged_pri_alt_ch }
    //
    // MODULE: MERGE PURGED CONTIGS AND HAPLOTIGS INTO ONE FILE
    //
    CAT_CAT_PURGEDUPS( purged_pri_alt_ch )

    //
    // LOGIC: DEFINE MERGED ASSEMBLY
    //
    merged_pri_alt = CAT_CAT_PURGEDUPS.out.file_out

    if ( polishing_on ) {
        //
        // MODULE: INDEX FASTA FOR THE MERGED PRIMARY CONTIGS AND HAPLOTIGS
        //
        SAMTOOLS_FAIDX_PURGEDUPS( CAT_CAT_PURGEDUPS.out.file_out, [[],[]] )
        ch_versions = ch_versions.mix(SAMTOOLS_FAIDX_PURGEDUPS.out.versions)

        //
        // LOGIC: CREATE AN ASSEMBLY CHANNEL FOR POLISHING
        //
        CAT_CAT_PURGEDUPS.out.file_out.join( SAMTOOLS_FAIDX_PURGEDUPS.out.fai )
                                  .set{ reference_ch }

        //
        // LOGIC: REFACTOR ILLUMINA CHANNEL TO PASS IT INTO THE POLISHING SUBWORKFLOW
        //
        PREPARE_INPUT.out.illumina_10X.map{ meta, reads, kmers -> [reads] }
                        .set{ illumina_10X_ch }
        
        //
        // SUBWORKFLOW: POLISH THE PRIMARY AND ALT
        //
        POLISHING(reference_ch, illumina_10X_ch, bed_chunks_polishing)
        ch_versions = ch_versions.mix(POLISHING.out.versions)
   
        //
        // LOGIC: UPDATE MERGED ASSEMBLY
        // 
        merged_pri_alt = POLISHING.out.fasta   

        //
        // MODULE: EXTRACT THE NAMES OF THE PRIMARY CONTIGS
        //
        KEEP_SEQNAMES_PRIMARY(PURGE_DUPS.out.pri)
        ch_versions = ch_versions.mix(KEEP_SEQNAMES_PRIMARY.out.versions)

        //
        // MODULE: SEPARATE POLISHED PRIMARY CONTIGS
        //
        SEQTK_SUBSEQ_PRIMARY(POLISHING.out.fasta, KEEP_SEQNAMES_PRIMARY.out.seqlist)
        ch_versions = ch_versions.mix(SEQTK_SUBSEQ_PRIMARY.out.versions)

        //
        // LOGIC: UPDATE THE PRIMARY CONTIGS CHANNEL WITH THE POLISHED 
        //        PRIMARY CONTIGS
        //
        POLISHING.out.fasta.map{ meta, f -> [id: meta.id] }
                            .combine(SEQTK_SUBSEQ_PRIMARY.out.sequences)
                            .set{ primary_contigs_ch }
        
        //
        // MODULE: EXTRACT THE NAMES OF THE HAPLOTIGS
        //
        KEEP_SEQNAMES_HAPLOTIGS(haplotigs_ch)

        //
        // MODULE: SEPARATE THE POLSIHED HAPLOTIGS
        //
        SEQTK_SUBSEQ_HAPLOTIGS(POLISHING.out.fasta, KEEP_SEQNAMES_HAPLOTIGS.out.seqlist)

        //
        // LOGIC: UPDATE THE HAPLOTIGS CHANNEL WITH THE POLISHED HAPLOTIGS
        //
        POLISHING.out.fasta.map{ meta, f -> [id: meta.id] }
                            .combine(SEQTK_SUBSEQ_HAPLOTIGS.out.sequences)
                            .set{ haplotigs_contigs_ch }

        //
        // LOGIC: COMBINE PRI AND ALT POLISHED CONTIGS INTO A CHANNEL
        //
        primary_contigs_ch.join(haplotigs_contigs_ch)
                        .map{ meta, pri, alt -> [[id:meta.id], pri, alt]}
                        .set{ polished_asm_stats_input_ch }

        //
        // SUBWORKFLOW: CALCULATE ASSEMBLY STATISTICS FOR THE POLISHED
        //              ASSEMBLY
        //
        GENOME_STATISTICS_POLISHED( polished_asm_stats_input_ch, 
                       PREPARE_INPUT.out.busco,
                       GENOMESCOPE_MODEL.out.hist,
                       GENOMESCOPE_MODEL.out.ktab
        )
    }

    //
    // LOGIC: CREATE A CHANNEL FOR THE PATHS TO HIC DATA
    //
    PREPARE_INPUT.out.hic.map{ meta, crams, motif -> [meta, crams] }
                         .set{ crams_ch }

    //
    // SUBWORKFLOW: MAP HIC DATA TO THE PRIMARY ASSEMBLY
    //
    HIC_MAPPING ( primary_contigs_ch,crams_ch )
    ch_versions = ch_versions.mix(HIC_MAPPING.out.versions)

    //
    // SUBWORKFLOW: SCAFFOLD THE PRIMARY ASSEMBLY
    //
    SCAFFOLDING( HIC_MAPPING.out.bed, primary_contigs_ch, cool_bin )
    ch_versions = ch_versions.mix(SCAFFOLDING.out.versions)

    //
    // LOGIC: CREATE A CHANNEL FOR THE FINAL ASSEMBLY REPRESENTED BY
    //        THE SCAFFOLDS AND HAPLOTIGS
    //
    SCAFFOLDING.out.fasta.combine(haplotigs_ch)
                        .map{meta_s, fasta_s, meta_h, fasta_h -> [ meta_h, fasta_s, fasta_h ]}
                        .set{ stats_input_ch }
    
    //
    // SUBWORKFLOW: CALCULATE ASSEMBLY STATISTICS FOR THE FINAL ASSEMBLY
    //
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
