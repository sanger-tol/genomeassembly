//
// Based on https://github.com/sanger-tol/readmapping/blob/9121ae2458f24e8fff6337a2f90e7ace1416b27f/subworkflows/local/align_short.nf
// from Sanger readmapping pipeline by @priyanka-surana
//
// Align short read (HiC and Illumina) data against the genome
//

include { SAMTOOLS_FASTQ } from '../../modules/nf-core/samtools/fastq/main'
include { BWAMEM2_MEM    } from '../../modules/nf-core/bwamem2/mem/main'
include { MARKDUP_STATS  } from '../../subworkflows/local/markdup_stats'
include { BWAMEM2_INDEX  } from '../../modules/nf-core/bwamem2/index/main'

workflow ALIGN_SHORT {
    take:
    fasta    // channel: [ val(meta), /path/to/fasta ]
    reads    // channel: [ val(meta), /path/to/datafile ]

    main:
    ch_versions = Channel.empty()

    // Convert from CRAM to FASTQ
    SAMTOOLS_FASTQ ( reads, true )
    ch_versions = ch_versions.mix(SAMTOOLS_FASTQ.out.versions.first())

    BWAMEM2_INDEX ( fasta )
    ch_versions = ch_versions.mix(BWAMEM2_INDEX.out.versions)

    // Align Fastq to Genome
    BWAMEM2_MEM ( SAMTOOLS_FASTQ.out.interleaved, BWAMEM2_INDEX.out.index, [] )
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
