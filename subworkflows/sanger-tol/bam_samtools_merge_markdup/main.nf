include { SAMTOOLS_FAIDX   } from '../../../modules/nf-core/samtools/faidx/main'
include { SAMTOOLS_MERGE   } from '../../../modules/nf-core/samtools/merge/main'
include { SAMTOOLS_MARKDUP } from '../../../modules/nf-core/samtools/markdup/main'

workflow BAM_SAMTOOLS_MERGE_MARKDUP {

    take:
    ch_bam              // channel: [ val(meta), [ bam ] ]
    ch_assemblies       // channel: [ val(meta), fasta ]
    val_mark_duplicates // boolean: mark duplicates on output bam

    main:
    ch_versions = Channel.empty()

    //
    // Module: Index assembly fastas
    //
    SAMTOOLS_FAIDX(
        ch_assemblies, // reference
        [ [:],[] ],    // fai
        false          // get sizes
    )
    ch_versions = ch_versions.mix(SAMTOOLS_FAIDX.out.versions)

    //
    // Logic: create a channel with both fai and gzi for each assembly
    //        We do it here so we don't cause downstream issues with the
    //        remainder join
    //
    ch_fai_gzi = SAMTOOLS_FAIDX.out.fai
        | join(SAMTOOLS_FAIDX.out.gzi, by: 0, remainder: true)
        | map { meta, fai, gzi -> [ meta, fai, gzi ?: [] ] }

    //
    // Logic: Prepare input for merging bams.
    //        We use the ch_chunk_counts to set a groupKey so that
    //        we emit groups downstream ASAP once all bams have been made
    //
    ch_samtools_merge_input = ch_bam
        | combine(ch_assemblies, by: 0)
        | combine(ch_fai_gzi, by: 0)
        | multiMap { meta, bams, assembly, fai, gzi ->
            bam:   [ meta, bams ]
            fasta: [ meta, assembly ]
            fai:   [ meta, fai ]
            gzi:   [ meta, gzi ]
        }

    //
    // Module: Merge position-sorted bam files
    //
    SAMTOOLS_MERGE(
        ch_samtools_merge_input.bam,
        ch_samtools_merge_input.fasta,
        ch_samtools_merge_input.fai,
        ch_samtools_merge_input.gzi,
    )
    ch_versions = ch_versions.mix(SAMTOOLS_MERGE.out.versions)

    //
    // Module: Mark duplicates on the merged bam
    //
    ch_samtools_markdup_input = SAMTOOLS_MERGE.out.bam
        | combine(ch_assemblies, by: 0)
        | filter { val_mark_duplicates }
        | multiMap { meta, bam, assembly ->
            bam:      [ meta, bam ]
            assembly: [ meta, assembly ]
        }

    SAMTOOLS_MARKDUP(
        ch_samtools_markdup_input.bam,
        ch_samtools_markdup_input.assembly
    )
    ch_versions = ch_versions.mix(SAMTOOLS_MARKDUP.out.versions)

    emit:
    bam      = val_mark_duplicates ? SAMTOOLS_MARKDUP.out.bam : SAMTOOLS_MERGE.out.bam // channel: [ val(meta), [ bam ] ]
    versions = ch_versions                                                             // channel: [ versions.yml ]
}
