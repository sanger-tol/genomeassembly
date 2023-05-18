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
include { LONGRANGER_MKREF } from '../../modules/local/longranger/mkref/main'
include { LONGRANGER_ALIGN } from '../../modules/local/longranger/align/main'

workflow POLISHING {
    take:
    fasta_in  //tuple meta, fasta, fai
    reads_10X // file
    bed_chunks_polishing    //val

    main:
    ch_versions = Channel.empty()

    //
    // Polishing step 1: map reads to the reference
    //
    fasta_in.map{ meta, fasta, fai -> [meta, fasta] }
            .set{ fasta_ch }
    LONGRANGER_MKREF(fasta_ch)
    ch_versions = ch_versions.mix(LONGRANGER_MKREF.out.versions)

    LONGRANGER_ALIGN( LONGRANGER_MKREF.out.folder, reads_10X )
    ch_versions = ch_versions.mix(LONGRANGER_ALIGN.out.versions)

    //
    // Polishing step 2: apply freebayes consensus based on longranger alignments
    //
    // Split genome into chunks
    fasta_in.map{ meta, fasta, fai -> [meta, fai] }
           .set{chunks_ch}
    BED_CHUNKS (chunks_ch, bed_chunks_polishing)
    ch_versions = ch_versions.mix(BED_CHUNKS.out.versions)
    intervals_structured = BED_CHUNKS.out.coords.toList().transpose()
    LONGRANGER_ALIGN.out.bam.join(LONGRANGER_ALIGN.out.bai)
                            .set{ bam_ch }
    intervals_freebayes = bam_ch.combine(intervals_structured)
     .map{ meta, bam, bai, bed -> [ [id: bed.getSimpleName()], bam, bai, [], [], bed] }

    // In case the average coverage from Longranger is provided use it for defining 
    // max coverage cut-off then scatter Freebayes over the genome chunks
    fasta = fasta_in.collect{it[1]}
    fai = fasta_in.collect{it[2]}
    LONGRANGER_ALIGN.out.csv.collect{it[1]}
                            .set{summary}
    LONGRANGER_COVERAGE(summary)
    ch_versions = ch_versions.mix(LONGRANGER_COVERAGE.out.versions)
    FREEBAYES(intervals_freebayes, fasta, fai, [], [], [], LONGRANGER_COVERAGE.out.cov)
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
    MERGE_FREEBAYES(BCFTOOLS_SORT.out.vcf.map{ meta, vcf -> [[id: 'merged'], vcf]}.groupTuple(), 
                    BCFTOOLS_SORT.out.vcf.map{ meta, vcf -> [meta, []] })
    ch_versions = ch_versions.mix(MERGE_FREEBAYES.out.versions)

    // Normalize variants and index normalized vcf
    MERGE_FREEBAYES.out.vcf.map{ meta, vcf -> [meta.id.toString(), vcf]}
        .join(MERGE_FREEBAYES.out.tbi.map{ meta, tbi -> [meta.id.toString(), tbi]})
        .combine(fasta_in)
        .map{ id_norm, vcf, tbi, meta, fasta, fai -> [meta, vcf, tbi] }
        .set{ input_norm }
    BCFTOOLS_NORM(input_norm, fasta)
    ch_versions = ch_versions.mix(BCFTOOLS_NORM.out.versions)
    BCFTOOLS_INDEX_NORM(BCFTOOLS_NORM.out.vcf)
    ch_versions = ch_versions.mix(BCFTOOLS_INDEX_NORM.out.versions)

    // Generate consensus fasta file    
    BCFTOOLS_NORM.out.vcf
        .join(BCFTOOLS_INDEX_NORM.out.tbi, by: [0], remainder: true)
        .join(fasta_ch, by: [0], remainder: true)
        .set{ch_merge}
    BCFTOOLS_CONSENSUS(ch_merge)
    ch_versions = ch_versions.mix(BCFTOOLS_CONSENSUS.out.versions)

    emit:
    fasta = BCFTOOLS_CONSENSUS.out.fasta
    versions = ch_versions.ifEmpty(null)
}
