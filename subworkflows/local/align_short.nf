//
// Based on https://github.com/sanger-tol/readmapping/blob/5a2657f5274a2c436b3896497a1322f2557bd9a8/subworkflows/local/align_short.nf
// from Sanger readmapping pipeline by @priyanka-surana
//
// Align short read (HiC and Illumina) data against the genome
//

include { SAMTOOLS_FASTQ } from '../../modules/nf-core/samtools/fastq/main'
include { BWAMEM2_MEM    } from '../../modules/nf-core/bwamem2/mem/main'
include { MARKDUP_STATS  } from '../../subworkflows/local/markdup_stats'

workflow ALIGN_SHORT {
    take:
    reads // channel: [ val(meta), [ datafile ] ]
    index // channel: [ val(meta), /path/to/bwamem2/ ]
    fasta // channel: /path/to/fasta

    main:
    ch_versions = Channel.empty()

    // Convert from CRAM to FASTQ
    SAMTOOLS_FASTQ ( reads )
    ch_versions = ch_versions.mix(SAMTOOLS_FASTQ.out.versions.first())

    // Align Fastq to Genome
    BWAMEM2_MEM ( SAMTOOLS_FASTQ.out.fastq, index, [] )
    ch_versions = ch_versions.mix(BWAMEM2_MEM.out.versions.first())

    // Merge, markdup, convert, and stats
    MARKDUP_STATS ( BWAMEM2_MEM.out.bam, fasta )
    ch_versions = ch_versions.mix(MARKDUP_STATS.out.versions)

    emit:
    bed = MARKDUP_STATS.out.bed
    cram = MARKDUP_STATS.out.cram
    crai = MARKDUP_STATS.out.crai
    stats = MARKDUP_STATS.out.stats
    idxstats = MARKDUP_STATS.out.idxstats
    flagstat = MARKDUP_STATS.out.flagstat

    versions = ch_versions
}
