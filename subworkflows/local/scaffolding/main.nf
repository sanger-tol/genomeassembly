include { BAMTOBED_SORT                              } from '../../../modules/local/bamtobedsort/main.nf'
include { SAMTOOLS_FAIDX as SAMTOOLS_FAIDX_CONTIGS   } from '../../../modules/nf-core/samtools/faidx/main.nf'
include { SAMTOOLS_FAIDX as SAMTOOLS_FAIDX_SCAFFOLDS } from '../../../modules/nf-core/samtools/faidx/main.nf'
include { YAHS                                       } from '../../../modules/nf-core/yahs/main'

// include { JUICER_PRE                          } from '../../../modules/local/juicer_pre.nf'
// include { JUICER_TOOLS_PRE                    } from '../../../modules/local/juicer_tools_pre.nf'
// include { PREPARE_PRETEXTMAP_INPUT            } from '../../../modules/local/prepare_pretext_input.nf'
// include { CHROM_SIZES                         } from '../../../modules/local/chrom_sizes.nf'
// include { COOLER_CLOAD                               } from '../../../modules/nf-core/cooler/cload/main.nf'
// include { COOLER_ZOOMIFY                             } from '../../../modules/nf-core/cooler/zoomify/main.nf'
// include { PRETEXTMAP                                 } from '../../../modules/nf-core/pretextmap/main.nf'
// include { PRETEXTSNAPSHOT                            } from '../../../modules/nf-core/pretextsnapshot/main'

workflow SCAFFOLDING {
    take:
    ch_fasta      // [meta, assembly]
    ch_map        // [meta, bam/bed]
    val_cool_bin  // val: cooler cload parameter

    main:
    ch_versions = Channel.empty()

    ch_map_split = ch_map
        | branch { meta, map ->
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
    ch_versions = ch_versions.mix(CONTIGS_FAIDX.out.versions)


    //
    // Module: scaffold contigs with YaHS
    //
    ch_yahs_input = ch_bed
        | combine(ch_fasta, by: 0)
        | combine(SAMTOOLS_FAIDX_CONTIGS.out.fai, by: 0)

    YAHS(ch_yahs_input)
    ch_versions = ch_versions.mix(YAHS.out.versions)

    // //
    // // Module: Index scaffolds
    // //
    // SCAFFOLDS_FAIDX(YAHS.out.scaffolds_fasta, [[],[]])
    // ch_versions = ch_versions.mix(SCAFFOLDS_FAIDX.out.versions)

    // //
    // // LOGIC: KEEP META
    // //
    // bed_in.map{ meta, bed -> meta}.set{ch_meta}

    // //
    // // LOGIC: PREPARE CONTACT PAIRS FOR COOLER
    // //
    // YAHS.out.binary.join(YAHS.out.scaffolds_agp)
    //                 .combine(scaf_ref)
    //                 .combine(scaf_ref_fai)
    //                 .map{meta, binary, agp, fa, fai -> [[id:meta.id, hap_id:hap_id], binary, agp, fai]}
    //                 .set{ch_merge}

    // //
    // // MODULE: PREPARE INPUT FOR COOLER
    // //
    // JUICER_PRE(ch_merge)
    // ch_versions = ch_versions.mix(JUICER_PRE.out.versions)

    // //
    // // LOGIC: BIN CONTACT PAIRS
    // //
    // JUICER_PRE.out.pairs.join(bed_in_hap)
    //                     .combine(Channel.of(cool_bin))
    //                     .set{ch_juicer}

    // //
    // // MODULE: GENERATE SCAFFOLD SIZES
    // //
    // CHROM_SIZES(SCAFFOLDS_FAIDX.out.fai)
    // ch_versions = ch_versions.mix(CHROM_SIZES.out.versions)

    // //
    // // MODULE: GENERATE A MULTI-RESOLUTION COOLER FILE BY COARSENING
    // //
    // COOLER_CLOAD(ch_juicer, CHROM_SIZES.out.chrom_sizes)
    // ch_versions = ch_versions.mix(COOLER_CLOAD.out.versions)

    // //
    // // LOGIC: REFACTOR CHANNEL FOR ZOOMIFY
    // //
    // COOLER_CLOAD.out.cool.map{ meta, cools, cool_bin-> [meta, cools]}
    //     .set{ch_cool}

    // //
    // //  MODULE: ZOOM COOL TO MCOOL
    // //
    // COOLER_ZOOMIFY(ch_cool)
    // ch_versions = ch_versions.mix(COOLER_ZOOMIFY.out.versions)

    // //
    // // LOGIC: EXTRACT INDEX FILE
    // //
    // SCAFFOLDS_FAIDX.out.fai.map{ meta, fai -> fai }.set{fai}

    // //
    // // MODULE: COMBINE SCAFFOLDS SIZES AND PAIRS FOR PRETEXT
    // //
    // PREPARE_PRETEXTMAP_INPUT(JUICER_PRE.out.pairs, fai)
    // ch_versions = ch_versions.mix(PREPARE_PRETEXTMAP_INPUT.out.versions)

    // //
    // // MODULE: GENERATE PRETEXT MAP FROM UPDATED PAIRS
    // //
    // PRETEXTMAP(PREPARE_PRETEXTMAP_INPUT.out.pairs, [])
    // ch_versions = ch_versions.mix(PRETEXTMAP.out.versions)

    // //
    // // MODULE: GENERATE PNG FROM STANDARD PRETEXT
    // //
    // PRETEXTSNAPSHOT(PRETEXTMAP.out.pretext)
    // ch_versions = ch_versions.mix(PRETEXTSNAPSHOT.out.versions)

    // //
    // // MODULE: GENERATE HIC MAP
    // //
    // JUICER_TOOLS_PRE(JUICER_PRE.out.pairs, CHROM_SIZES.out.chrom_sizes, 'yahs_scaffolds')
    // ch_versions = ch_versions.mix(JUICER_TOOLS_PRE.out.versions)

    emit:
    fasta = YAHS.out.scaffolds_fasta
    agp   = YAHS.out.agp
    versions = ch_versions
}
