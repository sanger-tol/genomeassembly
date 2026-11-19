include { BEDTOOLS_BAMTOBEDSORT                      } from '../../../modules/sanger-tol/bedtools/bamtobedsort/main'
include { RAGTAG_SCAFFOLD                            } from '../../../modules/nf-core/ragtag/scaffold/main'
include { SAMTOOLS_FAIDX as SAMTOOLS_FAIDX_CONTIGS   } from '../../../modules/nf-core/samtools/faidx/main'
include { SAMTOOLS_FAIDX as SAMTOOLS_FAIDX_SCAFFOLDS } from '../../../modules/nf-core/samtools/faidx/main'
include { YAHS                                       } from '../../../modules/nf-core/yahs/main'
include { YAHS_MAKEPAIRSFILE                         } from '../../../modules/sanger-tol/yahs/makepairsfile/main'

workflow FASTA_BAM_SCAFFOLDING_RAGTAG_YAHS {
    take:
    ch_fasta // [meta, fasta]
    ch_hic_bam // [meta, bam]
    ch_ragtag_reference // [meta, fasta]

    main:
    //
    // Module: Convert BAM to name-sorted BED
    //
    BEDTOOLS_BAMTOBEDSORT(ch_hic_bam)

    //
    // Module: Index input assemblies
    //
    SAMTOOLS_FAIDX_CONTIGS(
        ch_fasta.map { meta, assembly -> [meta, assembly, []] },
        false,
    )

    //
    // Module: scaffold additional haplotypes against primary haplotype
    //
    ch_ragtag_inputs = ch_ragtag_reference
        .combine(ch_fasta, by: 0)
        .multiMap { meta, ragtag_reference, contigs ->
            ragtag_reference: [meta, ragtag_reference]
            ragtag_contigs: [meta, contigs]
        }

    RAGTAG_SCAFFOLD(
        ch_ragtag_inputs.ragtag_contigs,
        ch_ragtag_inputs.ragtag_reference,
        [[], []],
        [[], [], []],
    )

    //
    // Module: run YaHS on the RagTag scaffolds
    //
    ch_yahs_ragtag_inputs = ch_fasta
        .combine(SAMTOOLS_FAIDX_CONTIGS.out.fai, by: 0)
        .combine(BEDTOOLS_BAMTOBEDSORT.out.sorted_bed, by: 0)
        .combine(RAGTAG_SCAFFOLD.out.corrected_agp, by: 0)

    YAHS(ch_yahs_ragtag_inputs)

    //
    // Module: Index output scaffolds
    //
    SAMTOOLS_FAIDX_SCAFFOLDS(
        YAHS.out.scaffolds_fasta.map { meta, assembly -> [meta, assembly, []] },
        true,
    )

    //
    // Module: Make pairs file to build contact maps with
    //
    ch_pairs_input = SAMTOOLS_FAIDX_SCAFFOLDS.out.fai
        .combine(YAHS.out.scaffolds_agp, by: 0)
        .combine(SAMTOOLS_FAIDX_CONTIGS.out.fai, by: 0)
        .combine(BEDTOOLS_BAMTOBEDSORT.out.sorted_bed, by: 0)

    YAHS_MAKEPAIRSFILE(ch_pairs_input)

    emit:
    scaffolds_fasta      = YAHS.out.scaffolds_fasta
    scaffolds_agp        = YAHS.out.scaffolds_agp
    scaffolds_pairs      = YAHS_MAKEPAIRSFILE.out.pairs
    scaffolds_chromsizes = SAMTOOLS_FAIDX_SCAFFOLDS.out.sizes
    ragtag_scaffolds     = RAGTAG_SCAFFOLD.out.corrected_assembly
    ragtag_agp           = RAGTAG_SCAFFOLD.out.corrected_agp
    ragtag_stats         = RAGTAG_SCAFFOLD.out.corrected_stats
    ragtag_paf           = RAGTAG_SCAFFOLD.out.ragtag_paf
    ragtag_delta         = RAGTAG_SCAFFOLD.out.ragtag_delta
    confidence_txt       = RAGTAG_SCAFFOLD.out.confidence_txt
    yahs_bin             = YAHS.out.binary
    yahs_inital          = YAHS.out.initial_break_agp
    yahs_intermediate    = YAHS.out.round_agp
    yahs_log             = YAHS.out.log
}
