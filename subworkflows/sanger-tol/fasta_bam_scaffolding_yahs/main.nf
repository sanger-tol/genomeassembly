include { BEDTOOLS_BAMTOBEDSORT                      } from '../../../modules/sanger-tol/bedtools/bamtobedsort/main'
include { SAMTOOLS_FAIDX as SAMTOOLS_FAIDX_CONTIGS   } from '../../../modules/nf-core/samtools/faidx/main'
include { SAMTOOLS_FAIDX as SAMTOOLS_FAIDX_SCAFFOLDS } from '../../../modules/nf-core/samtools/faidx/main'
include { YAHS                                       } from '../../../modules/nf-core/yahs/main'
include { YAHS_MAKEPAIRSFILE                         } from '../../../modules/sanger-tol/yahs/makepairsfile/main'

workflow FASTA_BAM_SCAFFOLDING_YAHS {
    take:
    ch_fasta // [meta, assembly]
    ch_hic_bam // [meta, bam/bed]

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
    // Module: scaffold contigs with YaHS
    //
    ch_yahs_input = ch_fasta
        .combine(SAMTOOLS_FAIDX_CONTIGS.out.fai, by: 0)
        .combine(BEDTOOLS_BAMTOBEDSORT.out.sorted_bed, by: 0)
        .map { meta, fasta, fai, bed ->
            // add empty AGP input
            [meta, fasta, fai, bed, []]
        }

    YAHS(ch_yahs_input)

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
    scaffolds_fasta   = YAHS.out.scaffolds_fasta
    scaffolds_agp     = YAHS.out.scaffolds_agp
    scaffolds_pairs   = YAHS_MAKEPAIRSFILE.out.pairs
    yahs_bin          = YAHS.out.binary
    yahs_inital       = YAHS.out.initial_break_agp
    yahs_intermediate = YAHS.out.round_agp
    yahs_log          = YAHS.out.log
}
