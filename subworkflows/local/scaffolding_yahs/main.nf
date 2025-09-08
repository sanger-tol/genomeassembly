include { BAMTOBED_SORT                              } from '../../../modules/local/bamtobedsort/main.nf'
include { COOLER_CLOAD                               } from '../../../modules/nf-core/cooler/cload/main.nf'
include { COOLER_ZOOMIFY                             } from '../../../modules/nf-core/cooler/zoomify/main.nf'
include { GAWK as GAWK_PROCESS_PAIRS_FILE            } from '../../../modules/nf-core/gawk/main.nf'
include { JUICERTOOLS_PRE                            } from '../../../modules/nf-core/juicertools/pre/main'
include { MAKE_PAIRS_FILE                            } from '../../../modules/local/make_pairs_file/main.nf'
include { PRETEXTMAP                                 } from '../../../modules/nf-core/pretextmap/main.nf'
include { PRETEXTSNAPSHOT                            } from '../../../modules/nf-core/pretextsnapshot/main.nf'
include { SAMTOOLS_FAIDX as SAMTOOLS_FAIDX_CONTIGS   } from '../../../modules/nf-core/samtools/faidx/main.nf'
include { SAMTOOLS_FAIDX as SAMTOOLS_FAIDX_SCAFFOLDS } from '../../../modules/nf-core/samtools/faidx/main.nf'
include { YAHS                                       } from '../../../modules/nf-core/yahs/main'

workflow SCAFFOLDING_YAHS {
    take:
    ch_fasta      // [meta, assembly]
    ch_map        // [meta, bam/bed]
    val_cool_bin  // val: cooler cload parameter

    main:
    ch_versions = Channel.empty()

    ch_map_split = ch_map
        | branch { _meta, map ->
            bam: map.getExtension() == "bam"
            bed: map.getExtension() == "bed"
        }

    //
    // Module: If map provided as BAM file - convert to name-sorted BED
    //
    BAMTOBED_SORT(
        ch_map_split.bam
    )
    ch_versions = ch_versions.mix(BAMTOBED_SORT.out.versions)

    ch_bed = ch_map_split.bed
        | mix(BAMTOBED_SORT.out.sorted_bed)

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
        | combine(SAMTOOLS_FAIDX_CONTIGS.out.fai, by: 0)
        | combine(ch_bed, by: 0)

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
        | join(YAHS.out.scaffolds_agp, by: 0)
        | join(SAMTOOLS_FAIDX_CONTIGS.out.fai, by: 0)
        | join(YAHS.out.binary, by: 0)

    MAKE_PAIRS_FILE(ch_pairs_input)
    ch_versions = ch_versions.mix(MAKE_PAIRS_FILE.out.versions)

    //
    // Module: Build PretextMap
    //
    PRETEXTMAP(
        MAKE_PAIRS_FILE.out.pairs, // Pairs file
        [[], [], []]
    )
    ch_versions = ch_versions.mix(PRETEXTMAP.out.versions)

    //
    // Module: Make a PNG of the PretextMap for fast viz
    //
    PRETEXTSNAPSHOT(PRETEXTMAP.out.pretext)
    ch_versions = ch_versions.mix(PRETEXTSNAPSHOT.out.versions)

    //
    // Module: Generate a multi-resolution cooler file by coarsening
    //
    ch_cooler_input = MAKE_PAIRS_FILE.out.pairs
        | map { meta, pairs -> [meta, pairs, []] }

    COOLER_CLOAD(
        ch_cooler_input,
        SAMTOOLS_FAIDX_SCAFFOLDS.out.sizes,
        "pairs",
        val_cool_bin
    )
    ch_versions = ch_versions.mix(COOLER_CLOAD.out.versions)

    //
    // Module: Zoom cool to mcool
    //
    COOLER_ZOOMIFY(COOLER_CLOAD.out.cool)
    ch_versions = ch_versions.mix(COOLER_ZOOMIFY.out.versions)


    //
    // Module: process .pairs file to remove the chromsize lines as juicer_pre
    // does not like them
    //
    GAWK_PROCESS_PAIRS_FILE(
        MAKE_PAIRS_FILE.out.pairs,
        "${projectDir}/bin/pairs_remove_chromsizes.awk",
        false
    )
    ch_versions = ch_versions.mix(GAWK_PROCESS_PAIRS_FILE.out.versions)

    //
    // Module: Generate juicer .hic map
    //
    ch_juicertools_pre_chrom_sizes = SAMTOOLS_FAIDX_SCAFFOLDS.out.sizes
        | map { meta, sizes ->
            [ meta, [], sizes ]
        }

    JUICERTOOLS_PRE(
        GAWK_PROCESS_PAIRS_FILE.out.output,
        ch_juicertools_pre_chrom_sizes
    )
    ch_versions = ch_versions.mix(JUICERTOOLS_PRE.out.versions)

    emit:
    assemblies  = YAHS.out.scaffolds_fasta
    agp         = YAHS.out.scaffolds_agp
    pretext     = PRETEXTMAP.out.pretext
    pretext_png = PRETEXTSNAPSHOT.out.image
    cool        = COOLER_ZOOMIFY.out.mcool
    hic         = JUICERTOOLS_PRE.out.hic
    versions    = ch_versions
}
