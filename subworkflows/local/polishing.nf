include { BCFTOOLS_VIEW                               } from '../../modules/nf-core/bcftools/view/main' 
include { BCFTOOLS_CONSENSUS                          } from '../../modules/nf-core/bcftools/consensus/main' 
include { BCFTOOLS_NORM                               } from '../../modules/nf-core/bcftools/norm/main' 
include { BCFTOOLS_CONCAT                             } from '../../modules/nf-core/bcftools/concat/main' 
include { BCFTOOLS_SORT                               } from '../../modules/nf-core/bcftools/sort/main'     
include { BCFTOOLS_INDEX as BCFTOOLS_INDEX_FB         } from '../../modules/nf-core/bcftools/index/main'     
include { BCFTOOLS_INDEX as BCFTOOLS_INDEX_NORM       } from '../../modules/nf-core/bcftools/index/main'     
include { GATK4_MERGEVCFS as MERGE_FREEBAYES          } from '../../modules/nf-core/gatk4/mergevcfs/main'
include { FREEBAYES                                   } from '../../modules/nf-core/freebayes/main' 
include { BED_CHUNKS                                  } from '../../modules/local/bed_chunks'
include { LONGRANGER_COVERAGE                         } from '../../modules/local/longranger_coverage'
include { LONGRANGER_MKREF                            } from '../../modules/local/longranger/mkref/main'
include { LONGRANGER_ALIGN                            } from '../../modules/local/longranger/align/main'

workflow POLISHING {
    take:
    fasta_in  //tuple meta, fasta, fai
    reads_10X // file
    bed_chunks_polishing    //val

    main:
    ch_versions = Channel.empty()

    //
    // LOGIC: SEPARATE ASSEMBLY FILE INTO CHANNEL
    //
    fasta_in.map{ meta, fasta, fai -> [meta, fasta] }
            .set{ fasta_ch }
    //
    // MODULE: GENERATE INDICES 
    //
    LONGRANGER_MKREF(fasta_ch)
    ch_versions = ch_versions.mix(LONGRANGER_MKREF.out.versions)

    //
    // MODULE: MAP 10X READS TO THE REFERENCE
    //
    LONGRANGER_ALIGN( LONGRANGER_MKREF.out.folder, reads_10X )
    ch_versions = ch_versions.mix(LONGRANGER_ALIGN.out.versions)

    //
    // LOGIC: SEPARATE INDEX FILE INTO CHANNEL
    //
    // Split genome into chunks
    fasta_in.map{ meta, fasta, fai -> [meta, fai] }
           .set{chunks_ch}

    //
    // MODULE: SPLIT ASSEMBLY INTO CHUNKS 
    //
    BED_CHUNKS (chunks_ch, bed_chunks_polishing)
    ch_versions = ch_versions.mix(BED_CHUNKS.out.versions)

    //
    // LOGIC: TRANSFORM CHUNKS CHANNEL INTO LIST OF INTERVALS
    //
    intervals_structured = BED_CHUNKS.out.coords.toList().transpose()

    //
    // LOGIC: JOIN READ MAPPING BAM WITH ITS INDEX
    //
    LONGRANGER_ALIGN.out.bam.join(LONGRANGER_ALIGN.out.bai)
                            .set{ bam_ch }
    
    //
    // LOGIC: CREATE DATA STRUCTURE FOR SCATTERING
    //
    intervals_freebayes = bam_ch.combine(intervals_structured)
     .map{ meta, bam, bai, bed -> [ [id: bed.getSimpleName()], bam, bai, [], [], bed] }

    //
    // LOGIC: SEPARATE ASSEMBLY AND ITS INDEX INTO CHANNELS
    //
    fasta = fasta_in.collect{it[1]}
    fai = fasta_in.collect{it[2]}

    //
    // LOGIC: EXTRACT ALIGNMENT SUMMARY FROM LONGRANGER RESULTS
    //
    LONGRANGER_ALIGN.out.csv.collect{it[1]}
                            .set{summary}

    //
    // MODULE: EXTRACT COVERAGE INFORMATION
    //
    LONGRANGER_COVERAGE(summary)
    ch_versions = ch_versions.mix(LONGRANGER_COVERAGE.out.versions)

    //
    // MODULE: SCATTER FREEBAYES OVER THE CHUNKS
    //
    FREEBAYES(intervals_freebayes, fasta, fai, [], [], [], LONGRANGER_COVERAGE.out.cov)
    ch_versions = ch_versions.mix(FREEBAYES.out.versions)

    //
    // MODULE: INDEX FREEBAYES OUTPUT
    //
    BCFTOOLS_INDEX_FB(FREEBAYES.out.vcf)
    ch_versions = ch_versions.mix(BCFTOOLS_INDEX_FB.out.versions)

    //
    // LOGIC: REFACTOR AND COMBINE VCF CHANNELS FOR FURTHER PROCESSING 
    //
    FREEBAYES.out.vcf.map{ meta, vcf -> [meta.id.toString(), vcf]} 
            .join(BCFTOOLS_INDEX_FB.out.tbi.map {meta, tbi -> [meta.id.toString(), tbi]})
            .map{ id, vcf, tbi -> [[ id: id.toString()+'_view'], vcf, tbi ]}
            .set{ input_view }

    //
    // MODULE: FILTER FREEBAYES RESULTS
    //
    BCFTOOLS_VIEW(input_view, [], [], [])
    ch_versions = ch_versions.mix(BCFTOOLS_VIEW.out.versions)

    //
    // LOGIC: REFACTOR CHANNEL TO AVOID NAME COLLISION
    //
    input_sort = BCFTOOLS_VIEW.out.vcf.map{ meta, vcf -> [ [id: meta.id.toString()+'_sorted'], vcf ]}

    //
    // MODULE: SORT FILTERED VCF
    //
    BCFTOOLS_SORT(input_sort)
    ch_versions = ch_versions.mix(BCFTOOLS_SORT.out.versions)
    
    //
    // LOGIC: SEPARATE META INTO CHANNEL
    //
    meta_ch = fasta_in.collect{it[0]}

    //
    // MODULE: MERGE FREEBAYES RESULTS ON CHUNKS
    //
    MERGE_FREEBAYES(BCFTOOLS_SORT.out.vcf.combine(fasta_in)
                        .map{ meta, vcf, meta_fin, fa, fai -> [[id: meta_fin.id], vcf]}.groupTuple(),
                    [ [id:'merged'], [] ] )
    ch_versions = ch_versions.mix(MERGE_FREEBAYES.out.versions)

    //
    // LOGIC: REFACTOR AND COMBINE CHANNELS FOR FURTHER PROCESSING
    //
    MERGE_FREEBAYES.out.vcf.map{ meta, vcf -> [meta.id.toString(), vcf]}
        .join(MERGE_FREEBAYES.out.tbi.map{ meta, tbi -> [meta.id.toString(), tbi]})
        .combine(fasta_in)
        .map{ id_norm, vcf, tbi, meta, fasta, fai -> [meta, vcf, tbi] }
        .set{ input_norm }
    
    //
    // LOGIC: CREATE CHANNEL FROM REFERENCE FILE AND META
    //
    fasta_in.map{ meta, fasta, fai -> [meta, fasta] }
            .set{ fasta_meta_ch }

    //
    // MODULE: LEFT-ALIGN AND NORMALIZE INDELS
    //
    BCFTOOLS_NORM(input_norm, fasta_meta_ch)
    ch_versions = ch_versions.mix(BCFTOOLS_NORM.out.versions)

    //
    // MODULE: INDEX NORMALIZED VARIANTS
    //
    BCFTOOLS_INDEX_NORM(BCFTOOLS_NORM.out.vcf)
    ch_versions = ch_versions.mix(BCFTOOLS_INDEX_NORM.out.versions)

    //
    // LOGIC: JOIN VCF CHANNEL WITH ITS INDEX AND FASTA REFERENCE CHANNELS
    //
    BCFTOOLS_NORM.out.vcf
        .join(BCFTOOLS_INDEX_NORM.out.tbi, by: [0], remainder: true)
        .join(fasta_ch, by: [0], remainder: true)
        .set{ ch_merge }

    //
    // MODULE: GENERATE CONSENSUS FASTA FILE
    //
    BCFTOOLS_CONSENSUS(ch_merge)
    ch_versions = ch_versions.mix(BCFTOOLS_CONSENSUS.out.versions)

    emit:
    fasta = BCFTOOLS_CONSENSUS.out.fasta
    versions = ch_versions.ifEmpty(null)
}
