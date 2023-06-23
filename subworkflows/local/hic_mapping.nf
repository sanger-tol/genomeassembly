//
// Based on https://github.com/sanger-tol/treeval/blob/f58fa556af60e427b6416a448b7bfc18bfe4a8c9/subworkflows/local/hic_mapping.nf and https://github.com/sanger-tol/readmapping/blob/9121ae2458f24e8fff6337a2f90e7ace1416b27f/subworkflows/local/align_short.nf
// from Sanger treeval and readmapping pipelines
//
// This subworkflow takes an input fasta sequence and csv style list of hico cram file to return
// merged alignment file and corresponding stats
// Align short read (HiC and Illumina) data against the genome
//

include { BWAMEM2_INDEX                             } from '../../modules/nf-core/bwamem2/index/main'
include { SAMTOOLS_FAIDX                            } from '../../modules/nf-core/samtools/faidx/main'
include { SAMTOOLS_MERGE                            } from '../../modules/nf-core/samtools/merge/main'
include { GENERATE_CRAM_CSV                         } from '../../modules/local/generate_cram_csv'
include { CRAM_FILTER_ALIGN_BWAMEM2_FIXMATE_SORT    } from '../../modules/local/cram_filter_align_bwamem2_fixmate_sort'
include { MARKDUP_STATS                             } from '../../subworkflows/local/markdup_stats'

workflow HIC_MAPPING {

    take:
    reference_tuple     // Channel [ val(meta), path(file) ]
    hic_reads_path      // Channel [ val(meta), path(directory) ]

    main:
    ch_versions = Channel.empty()

    //
    // MODULE: Indexing on reference output the folder of indexing files
    //
    BWAMEM2_INDEX (reference_tuple)
    ch_versions         = ch_versions.mix(BWAMEM2_INDEX.out.versions)
    
    //
    // LOGIC: make channel of hic reads as input for GENERATE_CRAM_CSV
    //
    reference_tuple
        .join( hic_reads_path )
        .map { meta, ref, hic_reads_path ->
                tuple([ id: meta.id, single_end: true], hic_reads_path) }
        .set { get_reads_input }

    //
    // MODULE: generate a cram csv file containing the required parametres for CRAM_FILTER_ALIGN_BWAMEM2_FIXMATE_SORT
    //
    GENERATE_CRAM_CSV ( get_reads_input )
    ch_versions         = ch_versions.mix(GENERATE_CRAM_CSV.out.versions)



    //
    // LOGIC: organise all parametres into a channel for CRAM_FILTER_ALIGN_BWAMEM2_FIXMATE_SORT
    //
    ch_filtering_input  = GENERATE_CRAM_CSV.out.csv
                            .splitCsv()
                            .combine (reference_tuple)
                            .combine (BWAMEM2_INDEX.out.index)
                            .map{ cram_id, cram_info, ref_id, ref_dir, bwa_id, bwa_path ->
                                  tuple([ 
                                        id: cram_id.id
                                        ], 
                                    file(cram_info[0]),
                                    cram_info[1],
                                    cram_info[2],
                                    cram_info[3],
                                    cram_info[4],
                                    cram_info[5],
                                    cram_info[6],
                                    bwa_path.toString() + '/' + ref_dir.toString().split('/')[-1])                          
                            }

    //
    // MODULE: parallel proccessing bwa-mem2 alignment by given interval of containers from cram files
    //
    CRAM_FILTER_ALIGN_BWAMEM2_FIXMATE_SORT ( ch_filtering_input  )
    ch_versions         = ch_versions.mix(CRAM_FILTER_ALIGN_BWAMEM2_FIXMATE_SORT.out.versions)

    //
    // LOGIC: PREPARING BAMS FOR MERGE
    //
    CRAM_FILTER_ALIGN_BWAMEM2_FIXMATE_SORT.out.mappedbam
        .map{ meta, file ->
            tuple( file )
        }
        .collect()
        .map { file ->
            tuple (
                [
                id: file[0].toString().split('/')[-1].split('_')[0]  // Change to sample_id
                ],
                file
            )
        }
        .set { collected_files_for_merge } 

    //
    // LOGIC: GENERATE INDEX OF REFERENCE
    //
    SAMTOOLS_FAIDX( reference_tuple )

    //
    // LOGIC: PREPARING MERGE INPUT
    //
    reference_tuple
        .combine( SAMTOOLS_FAIDX.out.fai )
        .multiMap { ref_meta, ref_fa, ref_idx_meta, ref_idx ->
            reference:  ref_fa
            ref_idx:  ref_idx
        }
        .set { ref_files }

    //
    // MODULE: MERGE POSITION SORTED BAM FILES AND MARK DUPLICATES
    //
    SAMTOOLS_MERGE ( collected_files_for_merge, ref_files.reference, ref_files.ref_idx )
    ch_versions         = ch_versions.mix ( SAMTOOLS_MERGE.out.versions.first() )

    // Merge, markdup, convert, and stats
    MARKDUP_STATS ( SAMTOOLS_MERGE.out.bam, reference_tuple )
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
