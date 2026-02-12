include { BCFTOOLS_VIEW                         } from '../../../modules/nf-core/bcftools/view'
include { BCFTOOLS_CONSENSUS                    } from '../../../modules/nf-core/bcftools/consensus'
include { BCFTOOLS_NORM                         } from '../../../modules/nf-core/bcftools/norm'
include { BCFTOOLS_SORT                         } from '../../../modules/nf-core/bcftools/sort'
include { BCFTOOLS_INDEX as BCFTOOLS_INDEX_FB   } from '../../../modules/nf-core/bcftools/index'
include { BCFTOOLS_INDEX as BCFTOOLS_INDEX_NORM } from '../../../modules/nf-core/bcftools/index'
include { BEDCHUNKS_CREATE                      } from '../../../modules/sanger-tol/bedchunks/create'
include { CAT_CAT as CONCATENATE_ASSEMBLIES     } from '../../../modules/nf-core/cat/cat'
include { FREEBAYES                             } from '../../../modules/nf-core/freebayes/main'
include { GATK4_MERGEVCFS                       } from '../../../modules/nf-core/gatk4/mergevcfs'
include { LONGRANGER_MKREF                      } from '../../../modules/sanger-tol/longranger/mkref'
include { LONGRANGER_ALIGN                      } from '../../../modules/sanger-tol/longranger/align'
include { SAMTOOLS_FAIDX                        } from '../../../modules/nf-core/samtools/faidx'
include { SEQKIT_GREP as SEPARATE_HAPLOTYPES    } from '../../../modules/nf-core/seqkit/grep'

workflow FASTA_10X_POLISHING_LONGRANGER_FREEBAYES {
    take:
    ch_assemblies          // [meta, [hap1, hap2, ... hapn]]
    ch_illumina_10x_reads  // [meta, reads]
    val_sequences_per_polishing_chunk // integer: number of sequences per polishing chunk

    main:
    //
    // Logic: rolling check of assembly meta objects to detect duplicates
    //
    def val_asm_meta_list = Collections.synchronizedSet(new HashSet())

    ch_assemblies
        .map { meta, _sample ->
            if (!val_asm_meta_list.add(meta)) {
                error("Error: Duplicate meta object found in `ch_assemblies` in FASTA_10X_POLISHING_LONGRANGER_FREEBAYES: ${meta}")
            }
            meta
        }

    //
    // Module: Concatenate haps together for polishing
    //
    CONCATENATE_ASSEMBLIES(ch_assemblies)

    //
    // Module: Index combined assembly FASTA
    //
    SAMTOOLS_FAIDX(
        CONCATENATE_ASSEMBLIES.out.file_out.map { meta, assembly -> [meta, assembly, []] },
        false)

    ch_assemblies_with_index = CONCATENATE_ASSEMBLIES.out.file_out.join(SAMTOOLS_FAIDX.out.fai)

    //
    // Module: Generate references
    //
    LONGRANGER_MKREF(CONCATENATE_ASSEMBLIES.out.file_out)

    //
    // Module: map 10x reads to the merged assemblies
    //
    ch_longranger_input = LONGRANGER_MKREF.out.reference
        .combine(ch_illumina_10x_reads, by: 0)
        .multiMap { meta, reference, reads ->
            reads: [meta, reads]
            reference: [meta, reference]
        }

    LONGRANGER_ALIGN(
        ch_longranger_input.reads,
        ch_longranger_input.reference,
    )

    //
    // Logic: Extract coverage information from Longranger summary and
    //        join to assembly
    //
    ch_longranger_coverage = LONGRANGER_ALIGN.out.csv
        .map { meta, summary ->
            def rows = summary.splitCsv(header: true, sep: ",")
            [meta, rows[0].mean_depth.toFloat().round().toLong()]
        }

    //
    // Module: split assembly into chunks
    //
    BEDCHUNKS_CREATE(
        SAMTOOLS_FAIDX.out.fai,
        val_sequences_per_polishing_chunk
    )

    //
    // Logic: Generate inputs for freebayes
    //
    ch_freebayes_input = ch_assemblies_with_index
        .combine(LONGRANGER_ALIGN.out.bam, by: 0)
        .combine(LONGRANGER_ALIGN.out.bai, by: 0)
        .combine(BEDCHUNKS_CREATE.out.bed, by: 0)
        .combine(ch_longranger_coverage, by: 0)
        .transpose(by: 5) // one entry per bed file
        .multiMap { meta, fasta, fai, bam, bai, bed, cov ->
            def chunk    = bed.name =~ /\.(\d+)\.bed$/
            def meta_new = meta + [longranger_cov: cov, chunk_id: chunk[0][1]]
            bam        : [meta_new, bam, bai, [], [], bed]
            fasta      : [meta_new, fasta]
            fai        : [meta_new, fai]
            samples    : [meta_new, []]
            populations: [meta_new, []]
            cnv        : [meta_new, []]
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

    //
    // Module: Index Freebayes output
    //
    BCFTOOLS_INDEX_FB(FREEBAYES.out.vcf)

    //
    // Logic: Refactor and combine VCF channels for further processing
    //
    ch_bcftools_view_input = FREEBAYES.out.vcf.combine(BCFTOOLS_INDEX_FB.out.tbi, by: 0)

    //
    // Module: Filer Freebayes results
    //
    BCFTOOLS_VIEW(ch_bcftools_view_input, [], [], [])

    //
    // Module: Sort filtered VCFs
    //
    BCFTOOLS_SORT(BCFTOOLS_VIEW.out.vcf)

    //
    // Module: Merge Freebayes results on chunks
    //
    ch_gatk4_mergevcf_input = BCFTOOLS_SORT.out.vcf
        .map { meta, vcf ->
            return [meta - meta.subMap(["longranger_cov", "chunk_id"]), vcf]
        }
        .groupTuple(by: 0)
        .map { meta, vcf -> [meta, vcf.sort { f -> f.getName() }] }

    GATK4_MERGEVCFS(
        ch_gatk4_mergevcf_input,
        [[:], []]
    )

    //
    // Module: Left-align and normalize indels
    //
    ch_bcftools_norm_input = ch_assemblies_with_index
        .combine(GATK4_MERGEVCFS.out.vcf, by: 0)
        .combine(GATK4_MERGEVCFS.out.tbi, by: 0)
        .multiMap { meta, fasta, _fai, vcf, tbi ->
            vcf  : [meta, vcf, tbi]
            fasta: [meta, fasta]
        }

    BCFTOOLS_NORM(
        ch_bcftools_norm_input.vcf,
        ch_bcftools_norm_input.fasta
    )

    //
    // Module: Index normalised variants
    //
    BCFTOOLS_INDEX_NORM(BCFTOOLS_NORM.out.vcf)

    //
    // Module: Generate consensus FASTA file
    //
    ch_bcftools_consensus_input = ch_assemblies_with_index
        .combine(BCFTOOLS_NORM.out.vcf      , by: 0)
        .combine(BCFTOOLS_INDEX_NORM.out.tbi, by: 0)
        .map { meta, fasta, _fai, vcf, tbi ->
            [meta, vcf, tbi, fasta, []]
        }

    BCFTOOLS_CONSENSUS(ch_bcftools_consensus_input)

    //
    // Module: split the haplotypes back out using seqkit grep + regex
    //
    ch_polished_assemblies_to_separate = BCFTOOLS_CONSENSUS.out.fasta
        .combine(ch_assemblies.map { meta, assemblies -> [meta, assemblies.size()] }, by: 0)
        .flatMap { meta, asm, size ->
            return (1..size).collect { count ->
                [meta + [_hap: "hap" + count, _nhaps: size], asm]
            }
        }

    SEPARATE_HAPLOTYPES(
        ch_polished_assemblies_to_separate,
        []
    )

    //
    // Logic: Take the split assemblies and re-orgnaise into [meta, [hap1, hap2, ... hapn]] format
    //
    ch_assemblies_polished_split = SEPARATE_HAPLOTYPES.out.filter
        .map { meta, asm ->
            def meta_new = meta - meta.subMap(["_hap", "_nhaps"])
            [ meta_new, meta._hap, asm ]
        }
        .groupTuple(by: 0)
        .map { meta, hap_labels, asms ->
            def hap_map = [hap_labels, asms].transpose().collectEntries()
            [ meta, hap_map.sort().values() ]
        }

    emit:
    polished_assemblies = ch_assemblies_polished_split
    assembly_consensus = BCFTOOLS_CONSENSUS.out.fasta
    longranger_bam = LONGRANGER_ALIGN.out.bam
    longranger_bai = LONGRANGER_ALIGN.out.bai
    longranger_summary = LONGRANGER_ALIGN.out.csv
    merged_vcf = GATK4_MERGEVCFS.out.vcf
    merged_vcf_tbi = GATK4_MERGEVCFS.out.tbi
}
