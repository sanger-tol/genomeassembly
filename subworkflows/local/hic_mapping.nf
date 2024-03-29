//
// Based on https://github.com/sanger-tol/treeval/blob/f58fa556af60e427b6416a448b7bfc18bfe4a8c9/subworkflows/local/hic_mapping.nf and https://github.com/sanger-tol/readmapping/blob/9121ae2458f24e8fff6337a2f90e7ace1416b27f/subworkflows/local/align_short.nf
// from Sanger treeval and readmapping pipelines
//
// This subworkflow takes an input fasta sequence and csv style list of hico cram file to return
// merged alignment file and corresponding stats
// Align short read (HiC and Illumina) data against the genome
//

include { BWAMEM2_INDEX                                       } from '../../modules/nf-core/bwamem2/index/main'
include { SAMTOOLS_FAIDX as SAMTOOLS_FAIDX_HIC_MAPPING        } from '../../modules/nf-core/samtools/faidx/main'
include { SAMTOOLS_MERGE as SAMTOOLS_MERGE_HIC_MAPPING        } from '../../modules/nf-core/samtools/merge/main'
include { SAMTOOLS_MARKDUP as SAMTOOLS_MARKDUP_HIC_MAPPING    } from '../../modules/nf-core/samtools/markdup/main'

include { BAMTOBED_SORT                                       } from '../../modules/local/bamtobed_sort'
include { GENERATE_CRAM_CSV                                   } from '../../modules/local/generate_cram_csv'
include { CRAM_FILTER_ALIGN_BWAMEM2_FIXMATE_SORT              } from '../../modules/local/cram_filter_align_bwamem2_fixmate_sort'
include { CONVERT_STATS                                       } from '../../subworkflows/local/convert_stats'

workflow HIC_MAPPING {

    take:
    reference_tuple     // Channel [ val(meta), path(file) ]
    hic_reads_path      // Channel [ val(meta), path(directory) ]

    main:
    ch_versions = Channel.empty()

    //
    // MODULE: INDEX REFERENCE FASTA
    //
    BWAMEM2_INDEX (reference_tuple)
    ch_versions = ch_versions.mix(BWAMEM2_INDEX.out.versions)
    
    //
    // LOGIC: JOIN HIC READS AND REFERENCE INTO A NEW CHANNEL
    //
    reference_tuple
        .join( hic_reads_path )
        .map { meta, ref, hic_reads_path ->
                tuple([ id: meta.id, single_end: true], hic_reads_path) }
        .set { get_reads_input }

    //
    // MODULE: GENERATE A CSV LISTING CRAM CHUNKS 
    //
    GENERATE_CRAM_CSV ( get_reads_input )
    ch_versions = ch_versions.mix(GENERATE_CRAM_CSV.out.versions)

    //
    // LOGIC: REFACTOR CHANNELS TO GET INPUT FOR CRAM_FILTER_ALIGN_BWAMEM2_FIXMATE_SORT
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
                                    cram_info[1], // crai path
                                    cram_info[2], // chunk starting position
                                    cram_info[3], // chunk end position
                                    cram_info[4], // basename
                                    cram_info[5], // the number of chunk
                                    cram_info[6], // rgline
                                    bwa_path.toString() + '/' + ref_dir.toString().split('/')[-1])                          
                            }

    //
    // MODULE: PERFORM READ MAPPING IN PARALLEL MANNER USING CRAM INTERVALS
    //
    CRAM_FILTER_ALIGN_BWAMEM2_FIXMATE_SORT ( ch_filtering_input  )
    ch_versions = ch_versions.mix(CRAM_FILTER_ALIGN_BWAMEM2_FIXMATE_SORT.out.versions)

    //
    // LOGIC: PREPARE BAMS FOR MERGE
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
    SAMTOOLS_FAIDX_HIC_MAPPING( reference_tuple, [[],[]] )

    //
    // LOGIC: PREPARE MERGE INPUT
    //
    reference_tuple
        .combine( SAMTOOLS_FAIDX_HIC_MAPPING.out.fai )
        .multiMap { ref_meta, ref_fa, ref_idx_meta, ref_idx ->
            reference_meta: [ref_meta, ref_fa]
            reference: ref_fa
            ref_idx:  [ref_idx_meta, ref_idx]
        }
        .set { ref_files }

    //
    // MODULE: MERGE POSITION SORTED BAM FILES 
    //
    SAMTOOLS_MERGE_HIC_MAPPING ( collected_files_for_merge, ref_files.reference_meta, ref_files.ref_idx )
    ch_versions = ch_versions.mix ( SAMTOOLS_MERGE_HIC_MAPPING.out.versions.first() )

    //
    // MODULE: MARK DUPLICATES ON THE MERGED BAM
    //
    SAMTOOLS_MARKDUP_HIC_MAPPING ( SAMTOOLS_MERGE_HIC_MAPPING.out.bam, ref_files.reference )
    ch_versions = ch_versions.mix ( SAMTOOLS_MARKDUP_HIC_MAPPING.out.versions )

    //
    // MODULE: FILTER OUT DUPLICATE READS, CONVERT BAM TO BED AND SORT BED FILE
    //
    BAMTOBED_SORT( SAMTOOLS_MARKDUP_HIC_MAPPING.out.bam )
    ch_versions = ch_versions.mix( BAMTOBED_SORT.out.versions )

    //
    // LOGIC: GENERATE INPUT FOR STATS SUBWORKFLOW
    //
    SAMTOOLS_MARKDUP_HIC_MAPPING.out.bam
                .map { meta, bam -> [ meta, bam, [] ] }
                .set { ch_stat }
    
    //
    // SUBWORKFLOW: PRODUCE READ MAPPING STATS
    //
    CONVERT_STATS ( ch_stat, ref_files.reference )
    ch_versions  = ch_versions.mix( CONVERT_STATS.out.versions )

    emit:
    bed = BAMTOBED_SORT.out.sorted_bed

    versions = ch_versions

}
