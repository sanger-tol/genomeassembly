//
// Merge and Markdup all alignments at specimen level
// Convert to CRAM and calculate statistics
//
include { SAMTOOLS_SORT } from '../../modules/nf-core/modules/samtools/sort/main'
include { SAMTOOLS_VIEW as SAMTOOLS_VIEW_MARKDUP } from '../../modules/nf-core/modules/samtools/view/main'
include { MARKDUPLICATE } from '../../subworkflows/local/markduplicate'
include { CONVERT_STATS } from '../../subworkflows/local/convert_stats'
include { BEDTOOLS_BAMTOBED } from '../../modules/sanger-tol/nf-core-modules/bedtools/bamtobed/main'
include { BEDTOOLS_SORT } from '../../modules/sanger-tol/nf-core-modules/bedtools/sort/main'

workflow MARKDUP_STATS {
    take:
    aln // channel: [ val(meta), [ sam ] ]
    fasta // channel: /path/to/fasta

    main:
    ch_versions = Channel.empty()

    // Sort SAM and convert to BAM
    SAMTOOLS_SORT ( aln )
    ch_versions = ch_versions.mix(SAMTOOLS_SORT.out.versions.first())

    // Collect all BWAMEM2 output by sample name
    SAMTOOLS_SORT.out.bam
//    .map { meta, bam ->
//    new_meta = meta.clone()
//    new_meta.id = new_meta.id.split('_')[0..-2].join('_')
//    [ [id: new_meta.id, datatype: new_meta.datatype] , bam ]
//    }
    .groupTuple(by: [0])
    .set { ch_bams }

    // Mark duplicates
    MARKDUPLICATE ( ch_bams )
    ch_versions = ch_versions.mix(MARKDUPLICATE.out.versions)

    // Convert merged BAM to CRAM and calculate indices and statistics
    ch_stat = MARKDUPLICATE.out.bam.map { meta, bam -> [ meta, bam, [] ] }
   
    SAMTOOLS_VIEW_MARKDUP ( ch_stat, fasta )
    ch_versions = ch_versions.mix(SAMTOOLS_VIEW_MARKDUP.out.versions)

    BEDTOOLS_BAMTOBED( SAMTOOLS_VIEW_MARKDUP.out.bam )
    ch_versions = ch_versions.mix(BEDTOOLS_BAMTOBED.out.versions)

    BEDTOOLS_SORT( BEDTOOLS_BAMTOBED.out.bed, "sorted.bed" )
    ch_versions = ch_versions.mix(BEDTOOLS_SORT.out.versions)

    CONVERT_STATS ( ch_stat, fasta )
    ch_versions = ch_versions.mix(CONVERT_STATS.out.versions)

    emit:
    bed = BEDTOOLS_SORT.out.sorted
    cram = CONVERT_STATS.out.cram
    crai = CONVERT_STATS.out.crai
    stats = CONVERT_STATS.out.stats
    idxstats = CONVERT_STATS.out.idxstats
    flagstat = CONVERT_STATS.out.flagstat

    versions = ch_versions
}
