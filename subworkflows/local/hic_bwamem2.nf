#!/usr/bin/env nextflow

// This subworkflow takes an input fasta sequence and csv style list of hic cram file to return
// alignment files including .mcool, pretext and .hic.
// Input - Assembled genomic fasta file, cram file directory
// Output - .mcool, .pretext, .hic

//
// MODULE IMPORT BLOCK
//
include { BWAMEM2_INDEX                                   } from '../../modules/nf-core/bwamem2/index/main'
include { CRAM_FILTER_ALIGN_BWAMEM2_FIXMATE_SORT          } from '../../modules/local/cram_filter_align_bwamem2_fixmate_sort'
include { SAMTOOLS_MERGE                                  } from '../../modules/nf-core/samtools/merge/main'

workflow HIC_BWAMEM2 {
    take:
    reference_tuple     // Channel: tuple [ val(meta), path( file )      ]
    csv_ch

    main:
    ch_versions         = Channel.empty()
    mappedbam_ch        = Channel.empty()

    BWAMEM2_INDEX (
        reference_tuple
        )
    ch_versions         = ch_versions.mix( BWAMEM2_INDEX.out.versions )

    csv_ch
        .splitCsv()
        .combine ( reference_tuple )
        .combine ( BWAMEM2_INDEX.out.index )
        .map{ cram_id, cram_info, ref_id, ref_dir, bwa_id, bwa_path ->
            tuple([
                    id: cram_id.id
                    ],
                file(cram_info[0]),
                cram_info[1], // crai path
                cram_info[2], // chunk starting position
                cram_info[3], // chunk end position
                cram_info[4], // basename
                cram_info[5], // the number of chunk
                cram_info[6], // rgline
                bwa_path.toString() + '/' + ref_dir.toString().split('/')[-1]
            )
    }
    .set { ch_filtering_input }

    //
    // MODULE: map hic reads by 10,000 container per time using bwamem2
    //
    CRAM_FILTER_ALIGN_BWAMEM2_FIXMATE_SORT (
        ch_filtering_input

    )
    ch_versions         = ch_versions.mix( CRAM_FILTER_ALIGN_BWAMEM2_FIXMATE_SORT.out.versions )
    mappedbam_ch        = CRAM_FILTER_ALIGN_BWAMEM2_FIXMATE_SORT.out.mappedbam

    //
    // LOGIC: PREPARING BAMS FOR MERGE
    //
    mappedbam_ch
        .map{ meta, file ->
            tuple( file )
        }
        .collect()
        .map { file ->
            tuple (
                [
                id: file[0].toString().split('/')[-1].split('_')[0]
                ],
                file
            )
        }
        .set { collected_files_for_merge }

    emit:
    mappedbams           = collected_files_for_merge
    versions            = ch_versions
}
