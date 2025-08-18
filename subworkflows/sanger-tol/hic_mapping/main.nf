include { BWAMEM2_INDEX              } from '../../../modules/nf-core/bwamem2/index/main'
include { HICCRAMALIGN_CHUNKS        } from '../../../modules/sanger-tol/hiccramalign/chunks'
include { HICCRAMALIGN_BWAMEM2ALIGN  } from '../../../modules/sanger-tol/hiccramalign/bwamem2align'
include { HICCRAMALIGN_MINIMAP2ALIGN } from '../../../modules/sanger-tol/hiccramalign/minimap2align'
include { MINIMAP2_INDEX             } from '../../../modules/nf-core/minimap2/index/main'
include { SAMTOOLS_FAIDX             } from '../../../modules/nf-core/samtools/faidx/main'
include { SAMTOOLS_INDEX             } from '../../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_MARKDUP           } from '../../../modules/nf-core/samtools/markdup/main'
include { SAMTOOLS_MERGE             } from '../../../modules/nf-core/samtools/merge/main'

workflow HIC_MAPPING {

    take:
    ch_assemblies        // Channel [meta, assembly]
    ch_hic_cram          // Channel [meta, cram] OR [meta, [cram1, cram2, ..., cram_n]]
    val_aligner          // string: [either "bwamem2" or "minimap2"]
    val_cram_chunk_size  // integer: Number of CRAM slices per chunk for mapping
    val_mark_duplicates  // boolean: Mark duplicates on the BAM?

    main:
    ch_versions = Channel.empty()

    //
    // Logic: check if CRAM files are accompanied by an index
    //        Get indexes, and index those that aren't
    //
    ch_hic_cram_raw = ch_hic_cram
        | transpose()
        | branch { meta, cram ->
            def cram_file = file(cram, checkIfExists: true)
            def index = cram + ".crai"
            have_index: file(index).exists()
                return [ meta, cram_file, file(index, checkIfExists: true) ]
            no_index: true
                return [ meta, cram_file ]
        }

    //
    // Module: Index CRAM files without indexes
    //
    SAMTOOLS_INDEX(ch_hic_cram_raw.no_index)
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions)

    ch_hic_cram_indexed = ch_hic_cram_raw.have_index
        | mix(
            ch_hic_cram_raw.no_index.join(SAMTOOLS_INDEX.out.crai)
        )

    //
    // Module: Process the cram index files to determine how many
    //         chunks to split into for mapping
    //
    HICCRAMALIGN_CHUNKS(
        ch_hic_cram_indexed,
        val_cram_chunk_size
    )
    ch_versions = ch_versions.mix(HICCRAMALIGN_CHUNKS.out.versions)

    //
    // Logic: Count the total number of cram chunks for downstream grouping
    //
    ch_n_cram_chunks = HICCRAMALIGN_CHUNKS.out.cram_slices
        | map { _meta, _cram, _crai, chunkn, _slices -> chunkn }
        | collect
        | map { chunkns -> chunkns.size() }

    //
    // Logic: Begin alignment - fork depending on specified aligner
    //
    if(val_aligner == "bwamem2") {
        //
        // Module: Create bwa-mem2 index for assembly
        //
        BWAMEM2_INDEX(ch_assemblies)
        ch_versions = ch_versions.mix(BWAMEM2_INDEX.out.versions)

        ch_assemblies_with_reference = ch_assemblies
            | combine(BWAMEM2_INDEX.out.index, by: 0)

        ch_cram_chunks = HICCRAMALIGN_CHUNKS.out.cram_slices
            | transpose()
            | combine(ch_assemblies_with_reference)
            | map { _meta, cram, crai, chunkn, slices, meta_assembly, index, assembly ->
                [ meta_assembly, cram, crai, chunkn, slices, index, assembly ]
            }

        HICCRAMALIGN_BWAMEM2ALIGN(ch_cram_chunks)
        ch_versions = ch_versions.mix(HICCRAMALIGN_BWAMEM2ALIGN.out.versions)

        ch_mapped_bams = HICCRAMALIGN_BWAMEM2ALIGN.out.bam
    } else if(val_aligner == "minimap2") {
        //
        // MODULE: generate minimap2 mmi file
        //
        MINIMAP2_INDEX(ch_assemblies)
        ch_versions = ch_versions.mix(MINIMAP2_INDEX.out.versions)

        ch_cram_chunks = HICCRAMALIGN_CHUNKS.out.cram_slices
            | transpose()
            | combine(MINIMAP2_INDEX.out.index)
            | map { _meta, cram, crai, chunkn, slices, meta_assembly, index ->
                [ meta_assembly, cram, crai, chunkn, slices, index ]
            }

        HICCRAMALIGN_MINIMAP2ALIGN(ch_cram_chunks)
        ch_versions = ch_versions.mix(HICCRAMALIGN_MINIMAP2ALIGN.out.versions)

        ch_mapped_bams = HICCRAMALIGN_MINIMAP2ALIGN.out.bam
    } else {
        log.error("Unsupported aligner: ${val_aligner}")
    }

    //
    // Logic: Index assembly fastas
    //
    SAMTOOLS_FAIDX(
        ch_assemblies, // reference
        [[:],[]],   // fai
        false       // get sizes
    )
    ch_versions = ch_versions.mix(SAMTOOLS_FAIDX.out.versions)

    //
    // Logic: Prepare input for merging bams.
    //        We use the ch_n_cram_chunks to set a groupKey so that
    //        we emit groups downstream ASAP once all bams have been made
    //
    ch_samtools_merge_input = ch_mapped_bams
        | combine(ch_n_cram_chunks)
        | map { meta, bam, n_chunks ->
            def key = groupKey(meta, n_chunks)
            [key, bam]
        }
        | groupTuple()
        | map { key, bam -> [key.target, bam] } // Get meta back out of groupKey
        | join(ch_assemblies, by: 0)
        | join(SAMTOOLS_FAIDX.out.fai, by: 0)
        | join(SAMTOOLS_FAIDX.out.gzi, by: 0, remainder: true)
        | multiMap { meta, bams, assembly, fai, gzi ->
            bam:   [meta, bams]
            fasta: [meta, assembly]
            fai:   [meta, fai]
            gzi:   [meta, gzi ?: []]
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
        | filter { val_mark_duplicates }
        | combine(ch_assemblies, by: 0)
        | multiMap { meta, bam, assembly ->
            bam:      [meta, bam]
            assembly: [meta, assembly]
        }

    SAMTOOLS_MARKDUP(
        ch_samtools_markdup_input.bam,
        ch_samtools_markdup_input.assembly
    )
    ch_versions = ch_versions.mix(SAMTOOLS_MARKDUP.out.versions)

    emit:
    bam      = val_mark_duplicates ? SAMTOOLS_MARKDUP.out.bam : SAMTOOLS_MERGE.out.bam
    versions = ch_versions
}
