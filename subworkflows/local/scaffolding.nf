include { COOLER_CLOAD                        } from '../../modules/nf-core/cooler/cload/main.nf'
include { COOLER_ZOOMIFY                      } from '../../modules/nf-core/cooler/zoomify/main.nf'
include { SAMTOOLS_FAIDX  as CONTIGS_FAIDX    } from '../../modules/nf-core/samtools/faidx/main.nf'
include { SAMTOOLS_FAIDX  as SCAFFOLDS_FAIDX  } from '../../modules/nf-core/samtools/faidx/main.nf'
include { YAHS                                } from '../../modules/nf-core/yahs/main'
include { PRETEXTMAP                          } from '../../modules/nf-core/pretextmap/main.nf'
include { PRETEXTSNAPSHOT                     } from '../../modules/nf-core/pretextsnapshot/main'

include { JUICER_PRE                          } from '../../modules/local/juicer_pre.nf'
include { JUICER_TOOLS_PRE                    } from '../../modules/local/juicer_tools_pre.nf'
include { PREPARE_PRETEXTMAP_INPUT            } from '../../modules/local/prepare_pretext_input.nf'
include { CHROM_SIZES                         } from '../../modules/local/chrom_sizes.nf'

workflow SCAFFOLDING {
    take:
    bed_in   // tuple(meta, bed)
    fasta_in // tuple(meta, fasta)
    cool_bin  // val: cooler cload parameter
    hap_id    // val: hap1/hap2/empty

    main:
    ch_versions = Channel.empty()

    //
    // MODULE: INDEX INPUT ASSEMBLY
    //
    CONTIGS_FAIDX( fasta_in, [[],[]] )
    ch_versions = ch_versions.mix(CONTIGS_FAIDX.out.versions)

    //
    // LOGIC: SEPARATE INPUT CHANNELS FOR YAHS
    //
    CONTIGS_FAIDX.out.fai.join( fasta_in )
                    .map{ meta, fai, fasta -> fasta }
                    .set{ scaf_ref }
    CONTIGS_FAIDX.out.fai.join( fasta_in )
                    .map{ meta, fai, fasta -> fai }
                    .set{ scaf_ref_fai }

    //
    // LOGIC: MIX IN THE HAPLOTYPE ID TO CONTROL THE OUTPUT SUFFIX
    //
    bed_in.map{ meta, bed -> [[id:meta.id, hap_id:hap_id],bed] }
        .set{ bed_in_hap }
    //
    // MODULE: PERFORM SCAAFFOLDING WITH YAHS
    //
    YAHS( bed_in_hap , scaf_ref, scaf_ref_fai )
    ch_versions = ch_versions.mix(YAHS.out.versions)

    //
    // MODULE: INDEX SCAFFOLDS
    //
    SCAFFOLDS_FAIDX(YAHS.out.scaffolds_fasta, [[],[]])
    ch_versions = ch_versions.mix(SCAFFOLDS_FAIDX.out.versions)

    //
    // LOGIC: KEEP META
    //
    bed_in.map{ meta, bed -> meta}.set{ch_meta}

    //
    // LOGIC: PREPARE CONTACT PAIRS FOR COOLER
    //
    YAHS.out.binary.join(YAHS.out.scaffolds_agp)
                    .combine(scaf_ref)
                    .combine(scaf_ref_fai)
                    .map{meta, binary, agp, fa, fai -> [[id:meta.id, hap_id:hap_id], binary, agp, fai]}
                    .set{ch_merge}

    //
    // MODULE: PREPARE INPUT FOR COOLER
    //
    JUICER_PRE(ch_merge)
    ch_versions = ch_versions.mix(JUICER_PRE.out.versions)

    //
    // LOGIC: BIN CONTACT PAIRS
    //
    JUICER_PRE.out.pairs.join(bed_in_hap)
                        .combine(Channel.of(cool_bin))
                        .set{ch_juicer}

    //
    // MODULE: GENERATE SCAFFOLD SIZES
    //
    CHROM_SIZES(SCAFFOLDS_FAIDX.out.fai)
    ch_versions = ch_versions.mix(CHROM_SIZES.out.versions)

    //
    // MODULE: GENERATE A MULTI-RESOLUTION COOLER FILE BY COARSENING
    //
    COOLER_CLOAD(ch_juicer, CHROM_SIZES.out.chrom_sizes)
    ch_versions = ch_versions.mix(COOLER_CLOAD.out.versions)

    //
    // LOGIC: REFACTOR CHANNEL FOR ZOOMIFY
    //
    COOLER_CLOAD.out.cool.map{ meta, cools, cool_bin-> [meta, cools]}
        .set{ch_cool}

    //
    //  MODULE: ZOOM COOL TO MCOOL
    //
    COOLER_ZOOMIFY(ch_cool)
    ch_versions = ch_versions.mix(COOLER_ZOOMIFY.out.versions)

    //
    // LOGIC: EXTRACT INDEX FILE
    //
    SCAFFOLDS_FAIDX.out.fai.map{ meta, fai -> fai }.set{fai}

    //
    // MODULE: COMBINE SCAFFOLDS SIZES AND PAIRS FOR PRETEXT
    //
    PREPARE_PRETEXTMAP_INPUT(JUICER_PRE.out.pairs, fai)
    ch_versions = ch_versions.mix(PREPARE_PRETEXTMAP_INPUT.out.versions)

    //
    // MODULE: GENERATE PRETEXT MAP FROM UPDATED PAIRS
    //
    PRETEXTMAP(PREPARE_PRETEXTMAP_INPUT.out.pairs, [])
    ch_versions = ch_versions.mix(PRETEXTMAP.out.versions)

    //
    // MODULE: GENERATE PNG FROM STANDARD PRETEXT
    //
    PRETEXTSNAPSHOT(PRETEXTMAP.out.pretext)
    ch_versions = ch_versions.mix(PRETEXTSNAPSHOT.out.versions)

    //
    // MODULE: GENERATE HIC MAP
    //
    JUICER_TOOLS_PRE(JUICER_PRE.out.pairs, CHROM_SIZES.out.chrom_sizes, 'yahs_scaffolds')
    ch_versions = ch_versions.mix(JUICER_TOOLS_PRE.out.versions)

    emit:
    fasta = YAHS.out.scaffolds_fasta

    versions = ch_versions
}
