//
// Based on https://github.com/sanger-tol/readmapping/blob/5a2657f5274a2c436b3896497a1322f2557bd9a8/subworkflows/local/convert_stats.nf
// from Sanger readmapping pipeline by @priyanka-surana
//

include { SAMTOOLS_INDEX_HIC_MAP } from '../../../modules/nf-core/samtools/index'
include { SAMTOOLS_STATS    } from '../../../modules/nf-core/samtools/stats'
include { SAMTOOLS_FLAGSTAT } from '../../../modules/nf-core/samtools/flagstat'
include { SAMTOOLS_IDXSTATS } from '../../../modules/nf-core/samtools/idxstats'

workflow HIC_MAPPING_STATS {
    take:
    bam         // [meta, bam]
    assemblies  // [meta, fasta]

    main:
    ch_versions = Channel.empty()

    //
    // Module: index bam
    //
    SAMTOOLS_INDEX(bam)
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions)

    ch_bam_bai = bam
        | combine(SAMTOOLS_INDEX.out.bai, by: 0)

    ch_stats_input = ch_bam_bai
        | combine(assemblies, by: 0)
        | multiMap { meta, in_bam, assembly ->
            bam:      [meta, in_bam]
            assembly: [meta, assembly]
        }

    //
    // Module: calculate stats
    //
    SAMTOOLS_STATS(ch_stats_input.bam, ch_stats_input.bam.assembly)
    ch_versions = ch_versions.mix(SAMTOOLS_STATS.out.versions)

    //
    // Module: Calculate statistics based on flag values
    //
    SAMTOOLS_FLAGSTAT(ch_bam_bai)
    ch_versions = ch_versions.mix(SAMTOOLS_FLAGSTAT.out.versions)

    //
    // Module: Calculate index statistics
    //
    SAMTOOLS_IDXSTATS(ch_bam_bai)
    ch_versions = ch_versions.mix(SAMTOOLS_IDXSTATS.out.versions)

    emit:
    versions = ch_versions
}
