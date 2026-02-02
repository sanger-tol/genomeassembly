include { BEDTOOLS_BAMTOBEDSORT                      } from '../../../modules/sanger-tol/bedtools/bamtobedsort/main'
include { SAMTOOLS_FAIDX as SAMTOOLS_FAIDX_CONTIGS   } from '../../../modules/nf-core/samtools/faidx/main'
include { SAMTOOLS_FAIDX as SAMTOOLS_FAIDX_SCAFFOLDS } from '../../../modules/nf-core/samtools/faidx/main'
include { YAHS                                       } from '../../../modules/nf-core/yahs/main'
include { YAHS_MAKEPAIRSFILE                         } from '../../../modules/sanger-tol/yahs/makepairsfile/main'

include { PAIRS_CREATE_CONTACT_MAPS                  } from '../pairs_create_contact_maps/main'

workflow FASTA_BAM_SCAFFOLDING_YAHS {

    take:
    ch_fasta          // [meta, assembly]
    ch_hic_bam        // [meta, bam/bed]
    val_build_pretext // bool: build pretext map
    val_build_cooler  // bool: build cooler
    val_build_juicer  // bool: build juicer
    val_cool_bin      // val: cooler cload parameter

    main:
    ch_versions = channel.empty()

    //
    // Module: Convert BAM to name-sorted BED
    //
    BEDTOOLS_BAMTOBEDSORT(ch_hic_bam)
    ch_versions = ch_versions.mix(BEDTOOLS_BAMTOBEDSORT.out.versions)

    //
    // Module: Index input assemblies
    //
    SAMTOOLS_FAIDX_CONTIGS(
        ch_fasta,
        [[],[]],
        false
    )
    ch_versions = ch_versions.mix(SAMTOOLS_FAIDX_CONTIGS.out.versions)

    //
    // Module: scaffold contigs with YaHS
    //
    ch_yahs_input = ch_fasta
        .combine(SAMTOOLS_FAIDX_CONTIGS.out.fai, by: 0)
        .combine(BEDTOOLS_BAMTOBEDSORT.out.sorted_bed, by: 0)
        .map { meta, fasta, fai, bed ->
            // add empty AGP input
            [ meta, fasta, fai, bed, [] ]
        }

    YAHS(ch_yahs_input)
    ch_versions = ch_versions.mix(YAHS.out.versions)

    //
    // Module: Index output scaffolds
    //
    SAMTOOLS_FAIDX_SCAFFOLDS(
        YAHS.out.scaffolds_fasta,
        [[],[]],
        true
    )
    ch_versions = ch_versions.mix(SAMTOOLS_FAIDX_SCAFFOLDS.out.versions)

    //
    // Module: Make pairs file to build contact maps with
    //
    ch_pairs_input = SAMTOOLS_FAIDX_SCAFFOLDS.out.fai
        .combine(YAHS.out.scaffolds_agp, by: 0)
        .combine(SAMTOOLS_FAIDX_CONTIGS.out.fai, by: 0)
        .combine(YAHS.out.binary, by: 0)

    YAHS_MAKEPAIRSFILE(ch_pairs_input)
    ch_versions = ch_versions.mix(YAHS_MAKEPAIRSFILE.out.versions)

    //
    // Subworkflow: Create Hi-C contact maps for visualisation of scaffolding outputs
    //
    ch_contact_map_inputs = YAHS_MAKEPAIRSFILE.out.pairs
        .combine(SAMTOOLS_FAIDX_SCAFFOLDS.out.sizes, by: 0)
        .multiMap { meta, pairs, sizes ->
            pairs: [ meta, pairs ]
            sizes: [ meta, sizes ]
        }

    PAIRS_CREATE_CONTACT_MAPS(
        ch_contact_map_inputs.pairs,
        ch_contact_map_inputs.sizes,
        val_build_pretext,
        val_build_cooler,
        val_build_juicer,
        val_cool_bin
    )
    ch_versions = ch_versions.mix(PAIRS_CREATE_CONTACT_MAPS.out.versions)

    emit:
    scaffolds_fasta   = YAHS.out.scaffolds_fasta
    scaffolds_agp     = YAHS.out.scaffolds_agp
    yahs_bin          = YAHS.out.binary
    yahs_inital       = YAHS.out.initial_break_agp
    yahs_intermediate = YAHS.out.round_agp
    yahs_log          = YAHS.out.log
    pretext           = PAIRS_CREATE_CONTACT_MAPS.out.pretext
    pretext_png       = PAIRS_CREATE_CONTACT_MAPS.out.pretext_png
    cool              = PAIRS_CREATE_CONTACT_MAPS.out.cool
    hic               = PAIRS_CREATE_CONTACT_MAPS.out.hic
    versions          = ch_versions
}
