/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//include { RAW_ASSEMBLY                                          } from '../subworkflows/local/raw_assembly'
//include { ORGANELLES                                            } from '../subworkflows/local/organelles'
include { KMERS                                     } from '../subworkflows/local/kmers'
//include { PURGE_DUPS                                            } from '../subworkflows/local/purge_dups'
//include { POLISHING                                             } from '../subworkflows/local/polishing'
//include { SCAFFOLDING                                           } from '../subworkflows/local/scaffolding'
//include { KEEP_SEQNAMES as KEEP_SEQNAMES_PRIMARY                } from '../modules/local/keep_seqnames'
//include { HIC_MAPPING                                           } from '../subworkflows/local/hic_mapping'
//include { GENOME_STATISTICS as GENOME_STATISTICS            } from '../subworkflows/local/genome_statistics'

//include { CAT_CAT as CAT_CAT_MITOHIFI_READS          } from "../modules/nf-core/cat/cat/main"
//include { SAMTOOLS_FAIDX as SAMTOOLS_FAIDX_PURGEDUPS } from '../modules/nf-core/samtools/faidx/main'
//include { SEQTK_SUBSEQ as SEQTK_SUBSEQ_HAPLOTIGS     } from '../modules/nf-core/seqtk/subseq/main'
//include { SEQTK_SUBSEQ as SEQTK_SUBSEQ_PRIMARY       } from '../modules/nf-core/seqtk/subseq/main'
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
    long_reads
    hic_reads
    illumina_10x
    mat_reads
    pat_reads
    busco
    mito
    plastid

    main:

    ch_versions = Channel.empty()

    //
    // SUBWORKFLOW: GENERATE KMER DATABASE AND PROFILE MODEL
    //
    KMERS(long_reads, mat_reads, pat_reads)
    ch_versions = ch_versions.mix(KMERS.out.versions)

//    if (params.hifiasm_trio_on) {
//
//        // LOGIC: PRODUCE YAK DATABASE FOR BOTH PAT AND MAT FOR HIFIASM
//        GENOMESCOPE_MODEL.out.patdb.map{ meta, patdb -> patdb }.set{ patdb_ch }
//        GENOMESCOPE_MODEL.out.matdb.map{ meta, matdb -> matdb }.set{ matdb_ch }
//
//        // LOGIC: PREPARE INPUT CHANNELS FOR GENOME_STATISTICS_SCAFFOLDS_TRIO
//        GENOMESCOPE_MODEL.out.phapktab
//        .map { it[1] }
//        .set { phapktab_ch }
//
//        GENOMESCOPE_MODEL.out.mhapktab
//        .map { it[1] }
//        .set { mhapktab_ch }
//
//        GENOMESCOPE_MODEL.out.phapktab
//        .map { it[2] }
//        .set { fastk_pktab}
//
//        GENOMESCOPE_MODEL.out.mhapktab
//        .map { it[2] }
//        .set { fastk_mktab}
//
//        // LOGIC: PREPARE INPUT CHANNELS FOR GENOME_STATISTICS_SCAFFOLDS_TRIO
//        RAW_ASSEMBLY( hifi_reads, hic_reads, params.hifiasm_hic_on, params.hifiasm_trio_on, patdb_ch, matdb_ch )
//    }
//    else {
//        RAW_ASSEMBLY( hifi_reads, hic_reads, params.hifiasm_hic_on, params.hifiasm_trio_on, [], [] )
//    }
//
//
//    //
//    // SUBWORKFLOW: RUN A HIFIASM ASSEMBLY ON THE HIFI READS; ALSO CREATE
//    //              A HIFIASM RUN IN HIC MODE IF THE FLAG IS SWITCHED ON
//    //
//    ch_versions = ch_versions.mix(RAW_ASSEMBLY.out.versions)
//
//    //
//    // LOGIC: DEFINE THE PRIMARY CONTIGS CHANNEL
//    //
//    RAW_ASSEMBLY.out.primary_contigs.set{ primary_contigs_ch }
//
//    //
//    // LOGIC: DEFINE THE HAPLOTIGS CHANNELS
//    //
//    RAW_ASSEMBLY.out.alternate_contigs.set{ haplotigs_ch }
//
//    //
//    // SUBWORKFLOW: CALCULATE STATISTICS FOR THE RAW ASSEMBLY
//    //
//    GENOME_STATISTICS_RAW( primary_contigs_ch.join(haplotigs_ch),
//        busco,
//        GENOMESCOPE_MODEL.out.hist,
//        GENOMESCOPE_MODEL.out.ktab,
//        [],
//        [],
//        [],
//        [],
//        false
//    )
//    ch_versions = ch_versions.mix(GENOME_STATISTICS_RAW.out.versions)
//
//    if ( params.organelles_on ) {
//        //
//        // LOGIC: CREATE CHANNEL FOR PRIMARY AND ALT CONTIGS
//        //
//        primary_contigs_ch.join(haplotigs_ch)
//            .map{ meta, pri, alt -> [meta, [pri, alt]]}
//            .set{ raw_pri_alt_ch }
//        //
//        // MODULE: MERGE PAW CONTIGS AND HAPLOTIGS INTO ONE FILE
//        //
//        CAT_CAT_RAW( raw_pri_alt_ch )
//
//        //
//        // LOGIC: DEFINE MERGED ASSEMBLY
//        //
//        merged_pri_alt_raw = CAT_CAT_RAW.out.file_out
//
//        //
//        // MODULE: MERGE INPUT FASTA FILES WITH PACBIO READS
//        //
//        CAT_CAT_MITOHIFI_READS(hifi_reads)
//        ch_versions = ch_versions.mix(CAT_CAT_MITOHIFI_READS.out.versions)
//
//        //
//        // SUBWORKFLOW: INDETIFY MITO IN THE RAW READS AND ASSEMBLY CONTIGS
//        //
//        ORGANELLES(CAT_CAT_MITOHIFI_READS.out.file_out, merged_pri_alt_raw,
//            mito, plastid)
//        ch_versions = ch_versions.mix(ORGANELLES.out.versions)
//    }
//
//    //
//    // LOGIC: CHECK IF THE HIFIASM HIC MODE WAS SWITCHED ON
//    //
//    if ( params.hifiasm_hic_on ) {
//        //
//        // SUBWORKFLOW: CALCULATE RAW ASSEMBLY STATISTICS FOR THE HIFIASN IN HIC MODE
//        //
//        GENOME_STATISTICS_RAW_HIC( RAW_ASSEMBLY.out.hap1_hic_contigs
//            .join(RAW_ASSEMBLY.out.hap2_hic_contigs),
//            busco,
//            GENOMESCOPE_MODEL.out.hist,
//            GENOMESCOPE_MODEL.out.ktab,
//            [],
//            [],
//            [],
//            [],
//            true
//        )
//    }
//
//    //
//    // LOGIC: CREATE AN INPUT DATA STRUCTURE FOR PURGING
//    //
//    hifi_reads.join(primary_contigs_ch)
//        .join(GENOMESCOPE_MODEL.out.model)
//        .set{ purge_dups_input }
//
//    //
//    // SUBWORKFLOW: RUN PURGE DUPS ON THE PRIMARY CONTIGS
//    //
//    if ( !params.hifiasm_trio_on ) {
//        PURGE_DUPS( purge_dups_input )
//        ch_versions = ch_versions.mix(PURGE_DUPS.out.versions)
//
//        //
//        // LOGIC: UPDATE THE PRIMARY CONTIGS CHANNEL
//        //
//        PURGE_DUPS.out.pri.map{ meta, fasta -> [[id:meta.id], fasta] }
//            .set{ primary_contigs_ch }
//
//        //
//        // LOGIC: SET APART THE HAPLOTIGS AFTER PURGING AND THE HIFIASM HAPLOTIGS
//        //
//        haplotigs_ch.combine( PURGE_DUPS.out.alt )
//            .map{ meta_h, h, meta_h_purged, h_purged -> [meta_h, [h, h_purged]]}
//            .set{ haplotigs_to_merge }
//
//        //
//        // MODULE: COMBINE PURGED SEQUENCES WITH THE ORIGINAL HAPLOTIGS
//        //
//        CAT_CAT_HAPLOTIGS{ haplotigs_to_merge }
//        haplotigs_ch = CAT_CAT_HAPLOTIGS.out.file_out
//
//        //
//        // SUBWORKFLOW: CALCULATE STATISTICS FOR THE PURGED ASSEMBLY
//        //
//        GENOME_STATISTICS_PURGED( primary_contigs_ch.join(haplotigs_ch),
//            busco,
//            GENOMESCOPE_MODEL.out.hist,
//            GENOMESCOPE_MODEL.out.ktab,
//            [],
//            [],
//            [],
//            [],
//            false
//        )
//
//        //
//        // LOGIC: CREATE A CHANNEL FOR THE PURGED CONTIGS AMD HAPLOTIGS
//        //
//        PURGE_DUPS.out.pri.join(haplotigs_ch)
//            .map{ meta, purged_pri, purged_alt -> [meta, [purged_pri, purged_alt]]}
//            .set{ purged_pri_alt_ch }
//        //
//        // MODULE: MERGE PURGED CONTIGS AND HAPLOTIGS INTO ONE FILE
//        //
//        CAT_CAT_PURGEDUPS( purged_pri_alt_ch )
//
//        //
//        // LOGIC: DEFINE MERGED ASSEMBLY
//        //
//        merged_pri_alt = CAT_CAT_PURGEDUPS.out.file_out
//    }
//
//    if ( params.polishing_on ) {
//        //
//        // MODULE: INDEX FASTA FOR THE MERGED PRIMARY CONTIGS AND HAPLOTIGS
//        //
//        SAMTOOLS_FAIDX_PURGEDUPS( CAT_CAT_PURGEDUPS.out.file_out, [[],[]] )
//        ch_versions = ch_versions.mix(SAMTOOLS_FAIDX_PURGEDUPS.out.versions)
//
//        //
//        // LOGIC: CREATE AN ASSEMBLY CHANNEL FOR POLISHING
//        //
//        CAT_CAT_PURGEDUPS.out.file_out.join( SAMTOOLS_FAIDX_PURGEDUPS.out.fai )
//            .set{ reference_ch }
//
//        //
//        // LOGIC: REFACTOR ILLUMINA CHANNEL TO PASS IT INTO THE POLISHING SUBWORKFLOW
//        //
//        illumina_10x.map{ meta, reads, kmers -> reads }
//            .set{ illumina_10X_ch }
//
//        //
//        // SUBWORKFLOW: POLISH THE PRIMARY AND ALT
//        //
//        POLISHING(reference_ch, illumina_10X_ch, params.bed_chunks_polishing)
//        ch_versions = ch_versions.mix(POLISHING.out.versions)
//
//        //
//        // LOGIC: UPDATE MERGED ASSEMBLY
//        //
//        merged_pri_alt = POLISHING.out.fasta
//
//        //
//        // MODULE: EXTRACT THE NAMES OF THE PRIMARY CONTIGS
//        //
//        KEEP_SEQNAMES_PRIMARY(PURGE_DUPS.out.pri)
//        ch_versions = ch_versions.mix(KEEP_SEQNAMES_PRIMARY.out.versions)
//
//        //
//        // MODULE: SEPARATE POLISHED PRIMARY CONTIGS
//        //
//        SEQTK_SUBSEQ_PRIMARY(POLISHING.out.fasta, KEEP_SEQNAMES_PRIMARY.out.seqlist)
//        ch_versions = ch_versions.mix(SEQTK_SUBSEQ_PRIMARY.out.versions)
//
//        //
//        // LOGIC: UPDATE THE PRIMARY CONTIGS CHANNEL WITH THE POLISHED
//        //        PRIMARY CONTIGS
//        //
//        POLISHING.out.fasta.map{ meta, f -> [id: meta.id] }
//            .combine(SEQTK_SUBSEQ_PRIMARY.out.sequences)
//            .set{ primary_contigs_ch }
//
//        //
//        // MODULE: EXTRACT THE NAMES OF THE HAPLOTIGS
//        //
//        KEEP_SEQNAMES_HAPLOTIGS(haplotigs_ch)
//
//        //
//        // MODULE: SEPARATE THE POLSIHED HAPLOTIGS
//        //
//        SEQTK_SUBSEQ_HAPLOTIGS(POLISHING.out.fasta, KEEP_SEQNAMES_HAPLOTIGS.out.seqlist)
//
//        //
//        // LOGIC: UPDATE THE HAPLOTIGS CHANNEL WITH THE POLISHED HAPLOTIGS
//        //
//        POLISHING.out.fasta.map{ meta, f -> [id: meta.id] }
//            .combine(SEQTK_SUBSEQ_HAPLOTIGS.out.sequences)
//            .set{ haplotigs_contigs_ch }
//
//        //
//        // LOGIC: COMBINE PRI AND ALT POLISHED CONTIGS INTO A CHANNEL
//        //
//        primary_contigs_ch.join(haplotigs_contigs_ch)
//            .map{ meta, pri, alt -> [[id:meta.id], pri, alt]}
//            .set{ polished_asm_stats_input_ch }
//
//        //
//        // SUBWORKFLOW: CALCULATE ASSEMBLY STATISTICS FOR THE POLISHED
//        //              ASSEMBLY
//        //
//        GENOME_STATISTICS_POLISHED( polished_asm_stats_input_ch,
//            busco,
//            GENOMESCOPE_MODEL.out.hist,
//            GENOMESCOPE_MODEL.out.ktab,
//            [],
//            [],  // Conditional pktab_ch
//            [],
//            [],
//            false
//        )
//    }
//
//    //
//    // LOGIC: CREATE A CHANNEL FOR THE PATHS TO HIC DATA
//    //
//    hic.map{ meta, crams, motif, hic_aligner -> [meta, crams] }
//        .set{ crams_ch }
//
//    hic.map{ meta, crams, motif, hic_aligner -> [meta, hic_aligner] }
//        .set{ hic_aligner_ch }
//
//    //
//    // SUBWORKFLOW: MAP HIC DATA TO THE PRIMARY ASSEMBLY
//    //
//    HIC_MAPPING ( primary_contigs_ch,crams_ch,hic_aligner_ch, "")
//    ch_versions = ch_versions.mix(HIC_MAPPING.out.versions)
//
//    //
//    // SUBWORKFLOW: SCAFFOLD THE PRIMARY ASSEMBLY
//    //
//    SCAFFOLDING( HIC_MAPPING.out.bed, primary_contigs_ch, params.cool_bin, "")
//    ch_versions = ch_versions.mix(SCAFFOLDING.out.versions)
//
//    //
//    // LOGIC: CREATE A CHANNEL FOR THE FINAL ASSEMBLY REPRESENTED BY
//    //        THE SCAFFOLDS AND HAPLOTIGS
//    //
//    SCAFFOLDING.out.fasta.combine(haplotigs_ch)
//                        .map{meta_s, fasta_s, meta_h, fasta_h -> [ meta_h, fasta_s, fasta_h ]}
//                        .set{ stats_input_ch }
//
//    //
//    // SUBWORKFLOW: CALCULATE ASSEMBLY STATISTICS FOR THE FINAL ASSEMBLY
//    //
//    GENOME_STATISTICS_SCAFFOLDS( stats_input_ch,
//        busco,
//        GENOMESCOPE_MODEL.out.hist,
//        GENOMESCOPE_MODEL.out.ktab,
//        [],
//        [],
//        [],
//        [],
//        false
//    )
//
//    if ( params.hifiasm_hic_on ) {
//        //
//        // SUBWORKFLOW: MAP HIC DATA TO THE HAP1 CONTIGS
//        //
//        HIC_MAPPING_HAP1 ( RAW_ASSEMBLY.out.hap1_hic_contigs, crams_ch, hic_aligner_ch, 'hap1' )
//        ch_versions = ch_versions.mix(HIC_MAPPING_HAP1.out.versions)
//
//        //
//        // SUBWORKFLOW: SCAFFOLD HAP1
//        //
//        SCAFFOLDING_HAP1( HIC_MAPPING_HAP1.out.bed, RAW_ASSEMBLY.out.hap1_hic_contigs, params.cool_bin, 'hap1' )
//        ch_versions = ch_versions.mix(SCAFFOLDING_HAP1.out.versions)
//
//        //
//        // SUBWORKFLOW: MAP HIC DATA TO THE HAP2 CONTIGS
//        //
//        HIC_MAPPING_HAP2 ( RAW_ASSEMBLY.out.hap2_hic_contigs, crams_ch, hic_aligner_ch, 'hap2' )
//        ch_versions = ch_versions.mix(HIC_MAPPING_HAP2.out.versions)
//
//        //
//        // SUBWORKFLOW: SCAFFOLD HAP2
//        //
//        SCAFFOLDING_HAP2( HIC_MAPPING_HAP2.out.bed, RAW_ASSEMBLY.out.hap2_hic_contigs, params.cool_bin, 'hap2' )
//        ch_versions = ch_versions.mix(SCAFFOLDING_HAP2.out.versions)
//
//        //
//        // LOGIC: CREATE A CHANNEL FOR THE FULL HAP1/HAP2 ASSEMBLY
//        //
//        SCAFFOLDING_HAP1.out.fasta.combine(SCAFFOLDING_HAP2.out.fasta)
//                    .map{meta_s, fasta_s, meta_h, fasta_h -> [ [id:meta_h.id], fasta_s, fasta_h ]}
//                    .set{ stats_haps_input_ch }
//
//        //
//        // SUBWORKFLOW: CALCULATE ASSEMBLY STATISTICS FOR HAP1/HAP2 ASSEMBLY
//        //
//        GENOME_STATISTICS_SCAFFOLDS_HAPS( stats_haps_input_ch,
//            busco,
//            GENOMESCOPE_MODEL.out.hist,
//            GENOMESCOPE_MODEL.out.ktab,
//            [],
//            [],
//            [],
//            [],
//            true
//        )
//    }
//
//    if ( params.hifiasm_trio_on ) {
//        //
//        // SUBWORKFLOW: MAP HIC DATA TO THE PATERNAL CONTIGS
//        //
//        HIC_MAPPING_PAT ( RAW_ASSEMBLY.out.pat_contigs, crams_ch, hic_aligner_ch, 'pat' )
//        ch_versions = ch_versions.mix(HIC_MAPPING_PAT.out.versions)
//
//        //
//        // SUBWORKFLOW: SCAFFOLD PAT
//        //
//        SCAFFOLDING_PAT( HIC_MAPPING_PAT.out.bed, RAW_ASSEMBLY.out.pat_contigs, params.cool_bin, 'pat' )
//        ch_versions = ch_versions.mix(SCAFFOLDING_PAT.out.versions)
//
//        //
//        // SUBWORKFLOW: MAP HIC DATA TO THE MATERNAL CONTIGS
//        //
//        HIC_MAPPING_MAT ( RAW_ASSEMBLY.out.mat_contigs, crams_ch, hic_aligner_ch, 'mat' )
//        ch_versions = ch_versions.mix(HIC_MAPPING_MAT.out.versions)
//
//        //
//        // SUBWORKFLOW: SCAFFOLD MAT
//        //
//        SCAFFOLDING_MAT( HIC_MAPPING_MAT.out.bed, RAW_ASSEMBLY.out.mat_contigs, params.cool_bin, 'mat' )
//        ch_versions = ch_versions.mix(SCAFFOLDING_MAT.out.versions)
//
//        //
//        // LOGIC: CREATE A CHANNEL FOR THE FULL PAT/MAT ASSEMBLY
//        //
//        SCAFFOLDING_PAT.out.fasta.combine(SCAFFOLDING_MAT.out.fasta)
//            .map{meta_s, fasta_s, meta_h, fasta_h -> [ [id:meta_h.id], fasta_s, fasta_h ]}
//            .set{ stats_trio_input_ch }
//
//        //
//        // SUBWORKFLOW: CALCULATE ASSEMBLY STATISTICS FOR TRIO ASSEMBLY
//        //
//        GENOME_STATISTICS_SCAFFOLDS_TRIO( stats_trio_input_ch,
//            busco,
//            GENOMESCOPE_MODEL.out.hist,
//            GENOMESCOPE_MODEL.out.ktab,
//            fastk_pktab,
//            phapktab_ch,
//            fastk_mktab,
//            mhapktab_ch,
//            true
//        )
//
//    }

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
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
