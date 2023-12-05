//
// Based on https://github.com/sanger-tol/readmapping/blob/5a2657f5274a2c436b3896497a1322f2557bd9a8/subworkflows/local/convert_stats.nf
// from Sanger readmapping pipeline by @priyanka-surana
//
// Convert BAM to CRAM, create index and calculate statistics
//

include { SAMTOOLS_VIEW     } from '../../modules/nf-core/samtools/view/main'
include { SAMTOOLS_INDEX    } from '../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_STATS    } from '../../modules/nf-core/samtools/stats/main'
include { SAMTOOLS_FLAGSTAT } from '../../modules/nf-core/samtools/flagstat/main'
include { SAMTOOLS_IDXSTATS } from '../../modules/nf-core/samtools/idxstats/main'

workflow CONVERT_STATS {
    take:
    bam
    fasta

    main:
    ch_versions = Channel.empty()

    // 
    // MODULE: CONVERT BAM TO CRAM
    //
    SAMTOOLS_VIEW ( bam, fasta, [] )
    ch_versions = ch_versions.mix(SAMTOOLS_VIEW.out.versions.first())

    //
    // MODULE: INDEX CRAM FILE
    //
    SAMTOOLS_INDEX ( SAMTOOLS_VIEW.out.cram )
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions.first())
    
    //
    // LOGIC: COMBINE CRAM AND CRAI INTP ONE CHANNEL
    //
    SAMTOOLS_VIEW.out.cram
    .join(SAMTOOLS_INDEX.out.crai, by: [0], remainder: true)
    .set { ch_cram_crai }

    //
    // MODULE: CALCULATE STATS
    //
    SAMTOOLS_STATS ( ch_cram_crai, fasta )
    ch_versions = ch_versions.mix(SAMTOOLS_STATS.out.versions.first())

    //
    // MODULE: CALCULATE STATISTTICS BASED ON FLAG VALUES
    //
    SAMTOOLS_FLAGSTAT ( ch_cram_crai )
    ch_versions = ch_versions.mix(SAMTOOLS_FLAGSTAT.out.versions.first())

    //
    // MODULE: CALCULATE INDEX STATISTICS
    //
    SAMTOOLS_IDXSTATS ( ch_cram_crai )
    ch_versions = ch_versions.mix(SAMTOOLS_IDXSTATS.out.versions.first())

    emit:
    versions = ch_versions
}
