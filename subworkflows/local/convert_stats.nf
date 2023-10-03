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

    // Convert BAM to CRAM
    SAMTOOLS_VIEW ( bam, fasta, [] )
    ch_versions = ch_versions.mix(SAMTOOLS_VIEW.out.versions.first())

    // Index CRAM file
    SAMTOOLS_INDEX ( SAMTOOLS_VIEW.out.cram )
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions.first())
    
    // Combine CRAM and CRAI into one channel
    SAMTOOLS_VIEW.out.cram
    .join(SAMTOOLS_INDEX.out.crai, by: [0], remainder: true)
    .set { ch_cram_crai }

    // Calculate statistics
    SAMTOOLS_STATS ( ch_cram_crai, fasta )
    ch_versions = ch_versions.mix(SAMTOOLS_STATS.out.versions.first())

    // Calculate statistics based on flag values
    SAMTOOLS_FLAGSTAT ( ch_cram_crai )
    ch_versions = ch_versions.mix(SAMTOOLS_FLAGSTAT.out.versions.first())

    // Calculate index statistics
    SAMTOOLS_IDXSTATS ( ch_cram_crai )
    ch_versions = ch_versions.mix(SAMTOOLS_IDXSTATS.out.versions.first())

    emit:
    cram = SAMTOOLS_VIEW.out.cram
    crai = SAMTOOLS_INDEX.out.crai
    stats = SAMTOOLS_STATS.out.stats
    flagstat = SAMTOOLS_FLAGSTAT.out.flagstat
    idxstats = SAMTOOLS_IDXSTATS.out.idxstats

    versions = ch_versions
}
