include { COOLER_CLOAD                               } from '../../../modules/nf-core/cooler/cload/main.nf'
include { COOLER_ZOOMIFY                             } from '../../../modules/nf-core/cooler/zoomify/main.nf'
include { GAWK as GAWK_PROCESS_PAIRS_FILE            } from '../../../modules/nf-core/gawk/main.nf'
include { JUICERTOOLS_PRE                            } from '../../../modules/nf-core/juicertools/pre/main'
include { PRETEXTMAP                                 } from '../../../modules/nf-core/pretextmap/main.nf'
include { PRETEXTSNAPSHOT                            } from '../../../modules/nf-core/pretextsnapshot/main.nf'

workflow PAIRS_CREATE_CONTACT_MAPS {
    take:
    ch_pairs          // [meta, pairs]
    ch_chrom_sizes    // [meta, sizes]
    val_build_pretext // bool: build pretext map
    val_build_cooler  // bool: build cooler
    val_build_juicer  // bool: build juicer
    val_cool_bin      // val: cooler cload parameter

    main:
    ch_versions = channel.empty()

    //
    // Module: Build PretextMap
    //
    PRETEXTMAP(
        ch_pairs.filter { val_build_pretext }, // Pairs file
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
    ch_cooler_input = ch_pairs
        .filter { val_build_cooler }
        .combine(ch_chrom_sizes, by: 0)
        .multiMap { meta, pairs, sizes ->
            pairs: [ meta, pairs, [] ]
            sizes: [ meta, sizes ]
        }

    COOLER_CLOAD(
        ch_cooler_input.pairs,
        ch_cooler_input.sizes,
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
    ch_pairs_remove_chromsizes_awk = channel.of('''\
        BEGIN { FS = OFS = "\\t" }
        !/^#chromsize/ {
            print $0
        }'''.stripIndent())
        .collectFile(name: "pairs_remove_chromsizes.awk", cache: true)
        .collect()

    GAWK_PROCESS_PAIRS_FILE(
        ch_pairs.filter { val_build_juicer },
        ch_pairs_remove_chromsizes_awk,
        false
    )
    ch_versions = ch_versions.mix(GAWK_PROCESS_PAIRS_FILE.out.versions)

    //
    // Module: Generate juicer .hic map
    //
    ch_juicertools_pre_input = GAWK_PROCESS_PAIRS_FILE.out.output
        .combine(ch_chrom_sizes, by: 0)
        .multiMap { meta, pairs, sizes ->
            pairs: [ meta, pairs ]
            sizes: [ meta, [], sizes ]
        }

    JUICERTOOLS_PRE(
        ch_juicertools_pre_input.pairs,
        ch_juicertools_pre_input.sizes
    )
    ch_versions = ch_versions.mix(JUICERTOOLS_PRE.out.versions)

    emit:
    pretext     = PRETEXTMAP.out.pretext
    pretext_png = PRETEXTSNAPSHOT.out.image
    cool        = COOLER_ZOOMIFY.out.mcool
    hic         = JUICERTOOLS_PRE.out.hic
    versions    = ch_versions
}
