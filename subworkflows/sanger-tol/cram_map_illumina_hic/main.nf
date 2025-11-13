include { BWAMEM2_INDEX              } from '../../../modules/nf-core/bwamem2/index/main'
include { CRAMALIGN_GENCRAMCHUNKS    } from '../../../modules/sanger-tol/cramalign/gencramchunks'
include { CRAMALIGN_BWAMEM2ALIGNHIC  } from '../../../modules/sanger-tol/cramalign/bwamem2alignhic'
include { CRAMALIGN_MINIMAP2ALIGNHIC } from '../../../modules/sanger-tol/cramalign/minimap2alignhic'
include { MINIMAP2_INDEX             } from '../../../modules/nf-core/minimap2/index/main'
include { SAMTOOLS_INDEX             } from '../../../modules/nf-core/samtools/index/main'

include { BAM_SAMTOOLS_MERGE_MARKDUP } from '../bam_samtools_merge_markdup/main'

workflow CRAM_MAP_ILLUMINA_HIC {

    take:
    ch_assemblies        // Channel [meta, assembly]
    ch_hic_cram          // Channel [meta, cram] OR [meta, [cram1, cram2, ..., cram_n]]
    val_aligner          // string: [either "bwamem2" or "minimap2"]
    val_cram_chunk_size  // integer: Number of CRAM slices per chunk for mapping

    main:
    ch_versions = channel.empty()

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
    ch_cram_indexes = ch_hic_cram_indexed
        | map { meta, cram, index -> [ meta, index ] }

    CRAMALIGN_GENCRAMCHUNKS(
        ch_cram_indexes,
        val_cram_chunk_size
    )
    ch_versions = ch_versions.mix(CRAMALIGN_GENCRAMCHUNKS.out.versions)

    //
    // Logic: Count the total number of cram chunks for downstream grouping
    //
    ch_n_cram_chunks = CRAMALIGN_GENCRAMCHUNKS.out.cram_slices
        | map { meta, chunkn, _slices -> [ meta, chunkn ] }
        | transpose()
        | groupTuple(by: 0)
        | map { meta, chunkns -> [ meta, chunkns.size() ] }

    //
    // Logic: Re-join the cram files and indexes to their chunk information
    //
    ch_cram_with_slices = ch_hic_cram_indexed
        | combine(CRAMALIGN_GENCRAMCHUNKS.out.cram_slices, by: 0)

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

        ch_cram_chunks = ch_cram_with_slices
            | transpose()
            | combine(ch_assemblies_with_reference, by: 0)

        CRAMALIGN_BWAMEM2ALIGNHIC(ch_cram_chunks)
        ch_versions = ch_versions.mix(CRAMALIGN_BWAMEM2ALIGNHIC.out.versions)

        ch_mapped_bams = CRAMALIGN_BWAMEM2ALIGNHIC.out.bam
    } else if(val_aligner == "minimap2") {
        //
        // MODULE: generate minimap2 mmi file
        //
        MINIMAP2_INDEX(ch_assemblies)
        ch_versions = ch_versions.mix(MINIMAP2_INDEX.out.versions)

        ch_cram_chunks = ch_cram_with_slices
            | transpose()
            | combine(MINIMAP2_INDEX.out.index, by: 0)

        CRAMALIGN_MINIMAP2ALIGNHIC(ch_cram_chunks)
        ch_versions = ch_versions.mix(CRAMALIGN_MINIMAP2ALIGNHIC.out.versions)

        ch_mapped_bams = CRAMALIGN_MINIMAP2ALIGNHIC.out.bam
    } else {
        log.error("Unsupported aligner: ${val_aligner}")
    }

    //
    // Logic: Prepare input for merging bams.
    //        We use the ch_n_cram_chunks to set a groupKey so that
    //        we emit groups downstream ASAP once all bams have been made
    //
    ch_merge_input = ch_mapped_bams
        | combine(ch_n_cram_chunks, by: 0)
        | map { meta, bam, n_chunks ->
            def key = groupKey(meta, n_chunks)
            [key, bam]
        }
        | groupTuple(by: 0)
        | map { key, bam -> [key.target, bam] } // Get meta back out of groupKey

    //
    // Subworkflow: merge BAM files and mark duplicates
    //
    BAM_SAMTOOLS_MERGE_MARKDUP(
        ch_merge_input,
        ch_assemblies,
        true
    )
    ch_versions = ch_versions.mix(BAM_SAMTOOLS_MERGE_MARKDUP.out.versions)

    emit:
    bam       = BAM_SAMTOOLS_MERGE_MARKDUP.out.bam
    bam_index = BAM_SAMTOOLS_MERGE_MARKDUP.out.bam_index
    versions  = ch_versions
}
