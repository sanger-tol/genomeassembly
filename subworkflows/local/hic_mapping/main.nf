//
// Based on https://github.com/sanger-tol/treeval/blob/f58fa556af60e427b6416a448b7bfc18bfe4a8c9/subworkflows/local/hic_mapping.nf and https://github.com/sanger-tol/readmapping/blob/9121ae2458f24e8fff6337a2f90e7ace1416b27f/subworkflows/local/align_short.nf
// from Sanger treeval and readmapping pipelines
//
// This subworkflow takes an input fasta sequence and csv style list of hico cram file to return
// merged alignment file and corresponding stats
// Align short read (HiC and Illumina) data against the genome
//

include { BAMTOBED_SORT                                    } from '../../../modules/local/bamtobed_sort'
include { BWAMEM2_INDEX                                    } from '../../../modules/nf-core/bwamem2/index/main'
include { HIC_MAPPING_STATS                                } from '../../../subworkflows/local/hic_mapping_stats'
include { CRAM_CHUNKS                                      } from '../../../modules/local/cramalign/chunks'
include { CRAM_FILTER_ALIGN_BWAMEM2_FIXMATE_SORT           } from '../../../modules/local/cramalign/cram_filter_align_bwamem2_fixmate_sort'
include { CRAM_FILTER_MINIMAP2_FILTER5END_FIXMATE_SORT     } from '../../../modules/local/cramalign/cram_filter_minimap2_filter5end_fixmate_sort'
include { MINIMAP2_INDEX                                   } from '../../../modules/nf-core/minimap2/index/main'
include { SAMTOOLS_FAIDX as SAMTOOLS_FAIDX_HIC_MAPPING     } from '../../../modules/nf-core/samtools/faidx/main'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_HIC_CRAM        } from '../../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_MARKDUP as SAMTOOLS_MARKDUP_HIC_MAPPING } from '../../../modules/nf-core/samtools/markdup/main'
include { SAMTOOLS_MERGE as SAMTOOLS_MERGE_HIC_MAPPING     } from '../../../modules/nf-core/samtools/merge/main'

workflow HIC_MAPPING {

    take:
    assemblies   // Channel [meta, assembly]
    hic_cram     // Channel [meta, crams]
    aligner
    hic_mapping_cram_chunk_size

    main:
    ch_versions = Channel.empty()

    //
    // Logic: check if CRAM files are accompanied by an index
    //        Get indexes, and index those that aren't
    //
    ch_hic_cram_raw = hic_cram
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
    SAMTOOLS_INDEX_HIC_CRAM(ch_hic_cram_raw.no_index)
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX_HIC_CRAM.out.versions)

    ch_hic_cram = ch_hic_cram_raw.have_index
        | mix(SAMTOOLS_INDEX_HIC_CRAM.out.crai)

    //
    // Module: Process the cram index files to determine how many
    //         chunks to split into for mapping
    //
    CRAM_CHUNKS(
        ch_hic_cram,
        hic_mapping_cram_chunk_size
    )

    //
    // Logic: Count the total number of cram chunks for downstream grouping
    //
    ch_n_cram_chunks = CRAM_CHUNKS.out.cram_slices
        | map { meta, cram, crai, chunkn, slices -> chunkn }
        | collect
        | map { chunkns -> chunkns.size() }

    //
    // Logic: Begin alignment - fork depending on specified aligner
    //
    if(aligner == "bwamem2") {
        //
        // Module: Create bwa-mem2 index for assembly
        //
        BWAMEM2_INDEX(assemblies)
        ch_versions = ch_versions.mix(BWAMEM2_INDEX.out.versions)

        ch_assemblies_with_reference = assemblies
            | combine(BWAMEM2_INDEX.out.index, by: 0)

        ch_cram_chunks = CRAM_CHUNKS.out.cram_slices
            | transpose()
            | combine(ch_assemblies_with_reference)
            | map { _meta, cram, crai, chunkn, slices, meta_assembly, index, assembly ->
                [meta_assembly, cram, crai, chunkn, slices, index, assembly]
            }

        CRAM_FILTER_ALIGN_BWAMEM2_FIXMATE_SORT(ch_cram_chunks)
        ch_versions = ch_versions.mix(CRAM_FILTER_ALIGN_BWAMEM2_FIXMATE_SORT.out.versions)

        ch_mapped_bams = CRAM_FILTER_ALIGN_BWAMEM2_FIXMATE_SORT.out.bam
    } else if (aligner == "minimap2") {
        //
        // MODULE: generate minimap2 mmi file
        //
        MINIMAP2_INDEX(assemblies)
        ch_versions = ch_versions.mix(MINIMAP2_INDEX.out.versions)

        ch_cram_chunks = CRAM_CHUNKS.out.cram_slices
            | transpose()
            | combine(MINIMAP2_INDEX.out.index)
            | map { _meta, cram, crai, chunkn, slices, meta_assembly, index ->
                [meta_assembly, cram, crai, chunkn, slices, index]
            }

        CRAM_FILTER_MINIMAP2_FILTER5END_FIXMATE_SORT(ch_cram_chunks)
        ch_versions = ch_versions.mix(CRAM_FILTER_MINIMAP2_FILTER5END_FIXMATE_SORT.out.versions)

        ch_mapped_bams = CRAM_FILTER_MINIMAP2_FILTER5END_FIXMATE_SORT.out.bam
    } else {
        log.error("Unsupported aligner: ${aligner}")
    }

    //
    // Logic: Index assembly fastas
    //
    SAMTOOLS_FAIDX_HIC_MAPPING(
        assemblies, // reference
        [[:],[]],   // fai
        false       // get sizes
    )

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
        | combine(assemblies, by: 0)
        | combine(SAMTOOLS_FAIDX_HIC_MAPPING.out.fai, by: 0)
        | multiMap { meta, bams, assembly, fai ->
            bam:   [meta, bams]
            fasta: [meta, fasta]
            fai:   [meta, fai]
        }

    //
    // Module: Merge position-sorted bam files
    //
    SAMTOOLS_MERGE_HIC_MAPPING(
        ch_samtools_merge_input.bam,
        ch_samtools_merge_input.fasta,
        ch_samtools_merge_input.fai
    )
    ch_versions = ch_versions.mix (SAMTOOLS_MERGE_HIC_MAPPING.out.versions)

    //
    // Module: Mark duplicates on the merged bam
    //
    ch_samtools_markdup_input = SAMTOOLS_MERGE_HIC_MAPPING.out.bam
        | combine(assemblies, by: 0)
        | multiMap { meta, bam, assembly ->
            bam:      [meta, bam]
            assembly: [meta, assembly]
        }

    SAMTOOLS_MARKDUP_HIC_MAPPING (
        ch_samtools_markdup_input.bam,
        ch_samtools_markdup_input.assembly
    )
    ch_versions = ch_versions.mix(SAMTOOLS_MARKDUP_HIC_MAPPING.out.versions)

    //
    // Module: Filter out duplicate reads convert bam to bed and sort BED file
    //
    BAMTOBED_SORT(SAMTOOLS_MARKDUP_HIC_MAPPING.out.bam)
    ch_versions = ch_versions.mix(BAMTOBED_SORT.out.versions)

    //
    // Subworkflow: Calculate stats for Hi-C mapping
    //
    HIC_MAPPING_STATS(
        SAMTOOLS_MARKDUP_HIC_MAPPING.out.bam,
        assemblies
    )
    ch_versions = ch_versions.mix(HIC_MAPPING_STATS.out.versions)

    emit:
    bam = SAMTOOLS_MARKDUP_HIC_MAPPING.out.bam
    bed = BAMTOBED_SORT.out.bed
}
