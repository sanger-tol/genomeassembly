//
// Based on https://github.com/sanger-tol/readmapping/blob/5a2657f5274a2c436b3896497a1322f2557bd9a8/subworkflows/local/markduplicate.nf
// from Sanger readmapping pipeline by @priyanka-surana
//
// Merge BAM files and mark duplicates
//
include { SAMTOOLS_MERGE    } from '../../modules/nf-core/samtools/merge/main'
include { SAMTOOLS_COLLATE  } from '../../modules/nf-core/samtools/collate/main'
include { SAMTOOLS_FIXMATE  } from '../../modules/nf-core/samtools/fixmate/main'
include { SAMTOOLS_SORT     } from '../../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_MARKDUP  } from '../../modules/nf-core/samtools/markdup/main'

workflow MARKDUPLICATE {
    take:
    bams // channel: [ val(meta), [ bams ] ]

    main:
    ch_versions = Channel.empty()

    // Merge position sorted bam files
    SAMTOOLS_MERGE ( bams, [], [] )
    ch_versions = ch_versions.mix(SAMTOOLS_MERGE.out.versions.first())

    // Collate merged BAM file
    SAMTOOLS_COLLATE ( SAMTOOLS_MERGE.out.bam, [] )
    ch_versions = ch_versions.mix(SAMTOOLS_COLLATE.out.versions.first())

    // Fill in mate coordinates and insert size fields
    SAMTOOLS_FIXMATE ( SAMTOOLS_COLLATE.out.bam )
    ch_versions = ch_versions.mix(SAMTOOLS_FIXMATE.out.versions.first())

    // Position sort BAM file
    SAMTOOLS_SORT ( SAMTOOLS_FIXMATE.out.bam )
    ch_versions = ch_versions.mix(SAMTOOLS_SORT.out.versions.first())

    // Mark duplicates
    SAMTOOLS_MARKDUP ( SAMTOOLS_SORT.out.bam, [] )
    ch_versions = ch_versions.mix(SAMTOOLS_MARKDUP.out.versions.first())

    emit:
    bam = SAMTOOLS_MARKDUP.out.bam

    versions = ch_versions
}
