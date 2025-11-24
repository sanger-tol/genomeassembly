include { SAMTOOLS_FAIDX    } from '../../../modules/nf-core/samtools/faidx/main'
include { SAMTOOLS_MERGE    } from '../../../modules/nf-core/samtools/merge/main'
include { SAMTOOLS_MERGEDUP } from '../../../modules/sanger-tol/samtools/mergedup/main'

workflow BAM_SAMTOOLS_MERGE_MARKDUP {

    take:
    ch_bam              // channel: [ val(meta), [ bam ] ]
    ch_assemblies       // channel: [ val(meta), fasta ]
    val_mark_duplicates // boolean: mark duplicates on output bam

    main:
    ch_versions = channel.empty()

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
    // Module: Either merge position-sorted bam files, or merge and mark duplicates
    //
    if(val_mark_duplicates) {
        SAMTOOLS_MERGEDUP(
            ch_samtools_merge_input.bam,
            ch_samtools_merge_input.fasta,
            ch_samtools_merge_input.fai,
            ch_samtools_merge_input.gzi,
        )
        ch_versions    = ch_versions.mix(SAMTOOLS_MERGEDUP.out.versions)

        ch_output_bam  = SAMTOOLS_MERGEDUP.out.bam
            | mix(SAMTOOLS_MERGEDUP.out.cram)

        ch_output_index = SAMTOOLS_MERGEDUP.out.csi
            | mix(SAMTOOLS_MERGEDUP.out.crai)

        ch_output_metrics = SAMTOOLS_MERGEDUP.out.metrics
    } else {
        SAMTOOLS_MERGE(
            ch_samtools_merge_input.bam,
            ch_samtools_merge_input.fasta,
            ch_samtools_merge_input.fai,
            ch_samtools_merge_input.gzi,
        )
        ch_versions    = ch_versions.mix(SAMTOOLS_MERGE.out.versions)

        ch_output_bam  = SAMTOOLS_MERGE.out.bam
            | mix(SAMTOOLS_MERGE.out.cram)

        ch_output_index = SAMTOOLS_MERGE.out.csi
            | mix(SAMTOOLS_MERGE.out.crai)

        ch_output_metrics = channel.empty()
    }

    emit:
    bam       = ch_output_bam     // channel: [ val(meta), path(bam) ]
    bam_index = ch_output_index   // channel: [ val(meta), path(index) ]
    metrics   = ch_output_metrics // channel [ val(meta), path(stats) ]
    versions  = ch_versions       // channel: [ versions.yml ]
}
