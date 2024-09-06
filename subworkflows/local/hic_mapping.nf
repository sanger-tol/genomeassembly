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

include { HIC_MINIMAP2                                        } from '../../subworkflows/local/hic_minimap2'
include { HIC_BWAMEM2                                         } from '../../subworkflows/local/hic_bwamem2'

workflow HIC_MAPPING {

    take:
    reference_tuple     // Channel [ val(meta), path(file) ]
    hic_reads_path      // Channel [ val(meta), path(directory) ]
    hic_aligner_ch      // Channel [ val(meta), val(hic_aligner)]
    hap_id              // Value hap_id

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
                tuple([ id: meta.id, hap_id: hap_id, single_end: true], hic_reads_path, hic_reads_path.collect { p -> p.resolveSibling(p.name + ".crai") } ) }
        .set { get_reads_input }

    //
    // MODULE: GENERATE A CSV LISTING CRAM CHUNKS 
    //
    GENERATE_CRAM_CSV ( get_reads_input )
    ch_versions = ch_versions.mix(GENERATE_CRAM_CSV.out.versions)

    hic_reads_path
        .combine(reference_tuple)
        .combine(hic_aligner_ch)
        .map{ meta, hic_read_path, ref_meta, ref, aligner_meta, hic_aligner ->
             tuple(
                [   id : ref_meta,
                    aligner : hic_aligner
                ],
                ref
             )
            }
        .branch{
            minimap2      : it[0].aligner == "minimap2"
            bwamem2       : it[0].aligner == "bwamem2"
        }
        .set{ch_aligner}

      
    //
    // SUBWORKFLOW: mapping hic reads using minimap2
    //
    HIC_MINIMAP2 (
        ch_aligner.minimap2,
        GENERATE_CRAM_CSV.out.csv
    )
    ch_versions         = ch_versions.mix( HIC_MINIMAP2.out.versions )
    mappedbams           = HIC_MINIMAP2.out.mappedbams
    
    //
    // SUBWORKFLOW: mapping hic reads using bwamem2
    //
    HIC_BWAMEM2 (
        ch_aligner.bwamem2,
        GENERATE_CRAM_CSV.out.csv
    )
    ch_versions         = ch_versions.mix( HIC_BWAMEM2.out.versions )
    mappedbams           = mappedbams.mix(HIC_BWAMEM2.out.mappedbams)

    mappedbams.map{meta, bams -> [[id: meta.id, hap_id:hap_id], bams]}
              .set { mappedbams }

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
    SAMTOOLS_MERGE_HIC_MAPPING ( mappedbams, ref_files.reference_meta, ref_files.ref_idx )
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
