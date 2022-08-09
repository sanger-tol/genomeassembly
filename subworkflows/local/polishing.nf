include { BCFTOOLS_VIEW                               } from '../../modules/nf-core/modules/bcftools/view/main' 
include { BCFTOOLS_CONSENSUS                          } from '../../modules/nf-core/modules/bcftools/consensus/main' 
include { BCFTOOLS_NORM                               } from '../../modules/nf-core/modules/bcftools/norm/main' 
include { BCFTOOLS_CONCAT                             } from '../../modules/nf-core/modules/bcftools/concat/main' 
include { BCFTOOLS_SORT                               } from '../../modules/nf-core/modules/bcftools/sort/main'     
include { BCFTOOLS_INDEX as BCFTOOLS_INDEX_FB         } from '../../modules/nf-core/modules/bcftools/index/main'     
include { BCFTOOLS_INDEX as BCFTOOLS_INDEX_NORM       } from '../../modules/nf-core/modules/bcftools/index/main'     
include { GATK4_MERGEVCFS as MERGE_FREEBAYES          } from '../../modules/nf-core/modules/gatk4/mergevcfs/main'
include { FREEBAYES                                   } from '../../modules/local/freebayes' 
include { BED_CHUNKS                                  } from '../../modules/local/bed_chunks'
include { LONGRANGER_COVERAGE                         } from '../../modules/local/longranger_coverage'

workflow POLISHING {
    take:
    bam_in //tuple meta, bam, bai
    fasta_in //tuple fasta, fai
    groups //val
    summary //path

    main:
    ch_versions = Channel.empty()
    meta = bam_in.collect{it[0]}
    fasta = fasta_in.collect{it[0]}
    fai = fasta_in.collect{it[1]}
    // Split genome into chunks
    bam_in.combine(fasta_in)
           .map{ meta, bam, bai, fasta, fai -> [meta, fai]}
           .set{chunks_ch}
    BED_CHUNKS (chunks_ch, groups)
    ch_versions = ch_versions.mix(BED_CHUNKS.out.versions)
    intervals_structured = BED_CHUNKS.out.coords.toList().transpose()
    intervals_freebayes = bam_in.combine(intervals_structured)
     .map{ meta, bam, bai, bed -> [ [id: bed.getSimpleName()], bam, bai, [], [], bed] }
    // In case the average coverage from Longranger is provided use it for defining 
    // max coverage cut-off then scatter Freebayes over the genome chunks
    if ( summary ) {
        LONGRANGER_COVERAGE(summary)
        ch_versions = ch_versions.mix(LONGRANGER_COVERAGE.out.versions)
        FREEBAYES(intervals_freebayes, fasta, fai, [], [], [], LONGRANGER_COVERAGE.out.cov)
    }
    else {
        FREEBAYES(intervals_freebayes, fasta, fai, [], [], [], [])
    }
    ch_versions = ch_versions.mix(FREEBAYES.out.versions)
    BCFTOOLS_INDEX_FB(FREEBAYES.out.vcf)
    ch_versions = ch_versions.mix(BCFTOOLS_INDEX_FB.out.versions)

    // Filter and sort vcf for each genome chunk
    FREEBAYES.out.vcf.map{ meta, vcf -> [meta.id.toString(), vcf]} 
            .join(BCFTOOLS_INDEX_FB.out.tbi.map {meta, tbi -> [meta.id.toString(), tbi]})
            .map{ id, vcf, tbi -> [[ id: id.toString()+'_view'], vcf, tbi ]}
            .set{ input_view }
    BCFTOOLS_VIEW(input_view, [], [], [])
    ch_versions = ch_versions.mix(BCFTOOLS_VIEW.out.versions)
    input_sort = BCFTOOLS_VIEW.out.vcf.map{ meta, vcf -> [ [id: meta.id.toString()+'_sorted'], vcf ]}
    BCFTOOLS_SORT(input_sort)
    ch_versions = ch_versions.mix(BCFTOOLS_SORT.out.versions)

    // Merge vcf files into one
    MERGE_FREEBAYES(BCFTOOLS_SORT.out.vcf.map{ meta, vcf -> [[id: 'merged'], vcf]}.groupTuple(), [])
    ch_versions = ch_versions.mix(MERGE_FREEBAYES.out.versions)

    // Normalize variants and index normalized vcf
    MERGE_FREEBAYES.out.vcf.map{ meta, vcf -> [meta.id.toString(), vcf]}
        .join(MERGE_FREEBAYES.out.tbi.map{ meta, tbi -> [meta.id.toString(), tbi]})
        .combine(bam_in)
        .map{ id_norm, vcf, tbi, meta, bam, bai -> [meta, vcf, tbi] }
        .set{ input_norm }
    BCFTOOLS_NORM(input_norm, fasta)
    ch_versions = ch_versions.mix(BCFTOOLS_NORM.out.versions)
    BCFTOOLS_INDEX_NORM(BCFTOOLS_NORM.out.vcf)
    ch_versions = ch_versions.mix(BCFTOOLS_INDEX_NORM.out.versions)

    // Generate consensus fasta file    
    fasta_in.combine(bam_in)
            .map{fasta, fai, meta, bam, bai -> [meta, fasta]}
            .set{ch_fasta}
    BCFTOOLS_NORM.out.vcf
        .join(BCFTOOLS_INDEX_NORM.out.tbi, by: [0], remainder: true)
        .join(ch_fasta, by: [0], remainder: true)
        .set{ch_merge}
    BCFTOOLS_CONSENSUS(ch_merge)
    ch_versions = ch_versions.mix(BCFTOOLS_CONSENSUS.out.versions)

    emit:
    fasta = BCFTOOLS_CONSENSUS.out.fasta
    versions = ch_versions.ifEmpty(null)
}
