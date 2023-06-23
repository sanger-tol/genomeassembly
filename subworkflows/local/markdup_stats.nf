//
// Based on https://github.com/sanger-tol/readmapping/blob/9121ae2458f24e8fff6337a2f90e7ace1416b27f/subworkflows/local/markdup_stats.nf
// from Sanger readmapping pipeline by @priyanka-surana
//
// Merge and Markdup all alignments at specimen level
// Convert to CRAM and calculate statistics
//
include { SAMTOOLS_SORT                           } from '../../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_VIEW as SAMTOOLS_VIEW_MARKDUP  } from '../../modules/nf-core/samtools/view/main'
include { MARKDUPLICATE                           } from '../../subworkflows/local/markduplicate'
include { CONVERT_STATS                           } from '../../subworkflows/local/convert_stats'
include { BEDTOOLS_BAMTOBED                       } from '../../modules/nf-core/bedtools/bamtobed/main'
include { GNU_SORT as BED_SORT                    } from '../../modules/local/gnu_sort'

workflow MARKDUP_STATS {
    take:
    aln // channel: [ val(meta), [ sam ] ]
    fasta // channel: [ val(meta), /path/to/fasta ]

    main:
    ch_versions = Channel.empty()

    // Sort SAM and convert to BAM
    SAMTOOLS_SORT ( aln )
    ch_versions = ch_versions.mix(SAMTOOLS_SORT.out.versions.first())

    // Collect all BWAMEM2 output by sample name
    SAMTOOLS_SORT.out.bam
    .groupTuple(by: [0])
    .set { ch_bams }

    // Mark duplicates
    MARKDUPLICATE ( ch_bams )
    ch_versions = ch_versions.mix(MARKDUPLICATE.out.versions)

    // Convert merged BAM to CRAM and calculate indices and statistics
    MARKDUPLICATE.out.bam
    | map { meta, bam -> [ meta, bam, [] ] }
    | set { ch_stat }

    fasta
    | map { meta, file -> file }
    | set { ch_fasta }
 
    SAMTOOLS_VIEW_MARKDUP ( ch_stat, ch_fasta, [] )
    ch_versions = ch_versions.mix(SAMTOOLS_VIEW_MARKDUP.out.versions)

    BEDTOOLS_BAMTOBED( SAMTOOLS_VIEW_MARKDUP.out.bam )
    ch_versions = ch_versions.mix(BEDTOOLS_BAMTOBED.out.versions)

    BED_SORT( BEDTOOLS_BAMTOBED.out.bed )
    ch_versions = ch_versions.mix(BED_SORT.out.versions)
    
    fasta
    | map { meta, file -> file }
    | set { ch_fasta }

    CONVERT_STATS ( ch_stat, ch_fasta )
    ch_versions = ch_versions.mix(CONVERT_STATS.out.versions)

    emit:
    bed = BED_SORT.out.bed
    cram = CONVERT_STATS.out.cram
    crai = CONVERT_STATS.out.crai
    stats = CONVERT_STATS.out.stats
    idxstats = CONVERT_STATS.out.idxstats
    flagstat = CONVERT_STATS.out.flagstat

    versions = ch_versions
}
