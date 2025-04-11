include { BCFTOOLS_VIEW                            } from '../../../modules/nf-core/bcftools/view'
include { BCFTOOLS_CONSENSUS                       } from '../../../modules/nf-core/bcftools/consensus'
include { BCFTOOLS_NORM                            } from '../../../modules/nf-core/bcftools/norm'
include { BCFTOOLS_CONCAT                          } from '../../../modules/nf-core/bcftools/concat'
include { BCFTOOLS_SORT                            } from '../../../modules/nf-core/bcftools/sort'
include { BCFTOOLS_INDEX as BCFTOOLS_INDEX_FB      } from '../../../modules/nf-core/bcftools/index'
include { BCFTOOLS_INDEX as BCFTOOLS_INDEX_NORM    } from '../../../modules/nf-core/bcftools/index'
include { CAT_CAT as CONCATENATE_ASSEMBLIES        } from '../../../modules/nf-core/cat/cat'
include { GATK4_MERGEVCFS as GATK4_MERGE_FREEBAYES } from '../../../modules/nf-core/gatk4/mergevcfs'
include { GAWK as GAWK_BED_CHUNKS                  } from '../../../modules/nf-core/gawk'
include { FREEBAYES                                } from '../../../modules/nf-core/freebayes/main'
include { LONGRANGER_MKREF                         } from '../../../modules/local/longranger/mkref'
include { LONGRANGER_ALIGN                         } from '../../../modules/local/longranger/align'
include { SAMTOOLS_FAIDX                           } from '../../../modules/nf-core/samtools/faidx'
include { SEQKIT_GREP as SEQKIT_GREP_SPLIT_HAPS    } from '../../../modules/nf-core/seqkit/grep/main'

workflow POLISHING_10X {
    take:
    assemblies             // [meta, fasta]
    illumina_10x_reads     // [meta, reads]

    main:
    ch_versions = Channel.empty()

    //
    // Module: concatenate hap1 and hap2 assemblies together
    //
    ch_assemblies_to_concatenate = assemblies
        | map { meta, asm ->
            def meta_new = meta - meta.subMap('haplotype')
            [meta_new, asm]
        }
        | groupTuple(by: 0, size: 2)

    CONCATENATE_ASSEMBLIES(ch_assemblies_to_concatenate)
    ch_versions = ch_versions.mix(CONCATENATE_ASSEMBLIES.out.versions)

    //
    // Module: Index merged assemblies
    //
    ch_merged_assemblies_to_index = CONCATENATE_ASSEMBLIES.out.file_out

    SAMTOOLS_FAIDX(
        ch_merged_assemblies_to_index,
        [[:], []],
        false)
    ch_versions = ch_versions.mix(SAMTOOLS_FAIDX.out.versions)

    ch_assemblies_with_index = SAMTOOLS_FAIDX.out.fa
        | join(SAMTOOLS_FAIDX.out.fai)

    //
    // Module: Generate references
    //
    LONGRANGER_MKREF(CONCATENATE_ASSEMBLIES.out.file_out)
    ch_versions = ch_versions.mix(LONGRANGER_MKREF.out.versions)

    //
    // Module: map 10x reads to the merged assemblies
    //
    LONGRANGER_ALIGN(
        LONGRANGER_MKREF.out.folder,
        illumina_10x_reads
    )
    ch_versions = ch_versions.mix(LONGRANGER_ALIGN.out.versions)

    //
    // Logic: Join read mapping BAM with its index
    //
    ch_aligned_bam = LONGRANGER_ALIGN.out.bam.join(LONGRANGER_ALIGN.out.bai)

    //
    // Logic: Extract coverage information from Longranger summary and
    //        join to assembly
    //
    ch_longranger_coverage = LONGRANGER_ALIGN.out.csv
        | map { meta, summary ->
            rows = summary.splitCsv(header: true, sep: ",")
            [meta, summary.mean_depth.toInteger()]
        }

    //
    // Module: split assembly into chunks
    //
    GAWK_BED_CHUNKS(SAMTOOLS_FAIDX.out.fai, [], true)
    ch_versions = ch_versions.mix(GAWK_BED_CHUNKS.out.versions)

    //
    // Logic: Generate inputs for freebayes
    //
    ch_freebayes_input = ch_assemblies_with_index
        | combine(LONGRANGER_ALIGN.out.bam  , by: 0)
        | combine(LONGRANGER_ALIGN.out.bai  , by: 0)
        | combine(GAWK_BED_CHUNKS.out.output, by: 0)
        | combine(ch_longranger_coverage, by: 0)
        | transpose(by: 5) // one entry per bed file
        | multiMap { meta, fasta, fai, bam, bai, bed, cov ->
            def chunk    = bed.name =~ /\.(\d+)\.bed$/
            def meta_new = meta + [longranger_cov: cov, chunk_id: chunk[0][1]]
            bam        : [meta_new, bam, bai, [], [], bed]
            fasta      : [meta_new, fasta                ]
            fai        : [meta_new, fai                  ]
            samples    : [meta_new, []                   ]
            populations: [meta_new, []                   ]
            cnv        : [meta_new, []                   ]
        }

    //
    // Module: Scatter Freebayes over the chunks
    //
    FREEBAYES(
        ch_freebayes_input.bam,
        ch_freebayes_input.fasta,
        ch_freebayes_input.fai,
        ch_freebayes_input.samples,
        ch_freebayes_input.populations,
        ch_freebayes_input.cnv
    )
    ch_versions = ch_versions.mix(FREEBAYES.out.versions)

    //
    // Module: Index Freebayes output
    //
    BCFTOOLS_INDEX_FB(FREEBAYES.out.vcf)
    ch_versions = ch_versions.mix(BCFTOOLS_INDEX_FB.out.versions)

    //
    // Logic: Refactor and combine VCF channels for further processing
    //
    ch_bcftools_view_input = FREEBAYES.out.vcf
        | combine(BCFTOOLS_INDEX_FB.out.tbi)

    //
    // MODULE: FILTER FREEBAYES RESULTS
    //
    BCFTOOLS_VIEW(ch_bcftools_view_input, [], [], [])
    ch_versions = ch_versions.mix(BCFTOOLS_VIEW.out.versions)

    //
    // Module: Sort filtered VCFs
    //
    BCFTOOLS_SORT(BCFTOOLS_VIEW.out.vcf)
    ch_versions = ch_versions.mix(BCFTOOLS_SORT.out.versions)

    //
    // Module: Merge Freebayes results on chunks
    //
    ch_merge_freebayes_input = BCFTOOLS_SORT.out.vcf
        | map { meta, vcf ->
            meta_new = meta - meta.subMap(["longranger_cov", "chunk_id"])
            [meta_new, vcf]
        }
        | groupTuple(by: 0)

    GATK4_MERGE_FREEBAYES(
        ch_merge_freebayes_input,
        [[:], []]
    )
    ch_versions = ch_versions.mix(GATK4_MERGE_FREEBAYES.out.versions)

    //
    // Module: Left-align and normalize indels
    //
    ch_bcftools_norm_input = ch_assemblies_with_index
        | combine(GATK4_MERGE_FREEBAYES.out.vcf, by: 0)
        | combine(GATK4_MERGE_FREEBAYES.out.tbi, by: 0)
        | multiMap{ meta, fasta, fai, vcf, tbi ->
            vcf  : [meta, vcf  ]
            fasta: [meta, fasta]
        }

    BCFTOOLS_NORM(
        ch_bcftools_norm_input.vcf,
        ch_bcftools_norm_input.fasta
    )
    ch_versions = ch_versions.mix(BCFTOOLS_NORM.out.versions)

    //
    // Module: Index normalised variants
    //
    BCFTOOLS_INDEX_NORM(BCFTOOLS_NORM.out.vcf)
    ch_versions = ch_versions.mix(BCFTOOLS_INDEX_NORM.out.versions)

    //
    // Module: Generate consensus FASTA file
    //
    ch_bcftools_consensus_input = ch_assemblies_with_index
        | combine(BCFTOOLS_NORM.out.vcf      , by: 0)
        | combine(BCFTOOLS_INDEX_NORM.out.tbi, by: 0)
        | map { meta, fasta, fai, vcf, tbi ->
            [meta, vcf, tbi, fasta, []]
        }

    BCFTOOLS_CONSENSUS(ch_bcftools_consensus_input)
    ch_versions = ch_versions.mix(BCFTOOLS_CONSENSUS.out.versions)

    //
    // Module: Separate back out primary/alt/hap1/hap2 contigs
    //
    ch_haps = Channel.of(["hap1", "hap2"])
    ch_assemblies_to_separate = BCFTOOLS_CONSENSUS.out.fasta
        | combine(ch_haps)
        | transpose(by: 2)
        | map { meta, asm, hap ->
            [meta + [haplotype: hap], asm]
        }

    SEQKIT_GREP_SPLIT_HAPS(
        ch_assemblies_to_separate,
        []
    )
    ch_versions         = ch_versions.mix(SEQKIT_GREP_SPLIT_HAPS.out.versions)
    ch_assemblies_split = SEQKIT_GREP_SPLIT_HAPS.out.filter
        | map { meta, fasta ->
            [meta + [assembly_stage: "polished"], fasta]
        }

    emit:
    assemblies = ch_assemblies_split
    versions   = ch_versions
}
