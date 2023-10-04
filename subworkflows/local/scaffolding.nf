include { COOLER_CLOAD                        } from '../../modules/nf-core/cooler/cload/main.nf' 
include { COOLER_ZOOMIFY                      } from '../../modules/nf-core/cooler/zoomify/main.nf' 
include { SAMTOOLS_FAIDX  as CONTIGS_FAIDX    } from '../../modules/nf-core/samtools/faidx/main.nf'
include { SAMTOOLS_FAIDX  as SCAFFOLDS_FAIDX  } from '../../modules/nf-core/samtools/faidx/main.nf'
include { YAHS                                } from '../../modules/nf-core/yahs/main'
include { JUICER_PRE                          } from '../../modules/local/juicer_pre.nf'
include { JUICER_TOOLS_PRE                    } from '../../modules/local/juicer_tools_pre.nf'
include { PREPARE_PRETEXTMAP_INPUT            } from '../../modules/local/prepare_pretext_input.nf'
include { PRETEXTMAP                          } from '../../modules/nf-core/pretextmap/main.nf'
include { PRETEXTSNAPSHOT                     } from '../../modules/nf-core/pretextsnapshot/main'
include { CHROM_SIZES                         } from '../../modules/local/chrom_sizes.nf'

workflow SCAFFOLDING {
    take:  
    bed_in   // tuple(meta, bed)
    fasta_in // tuple(meta, fasta)
    cool_bin  // val: cooler cload parameter 
    
    main:
    ch_versions = Channel.empty()
    CONTIGS_FAIDX( fasta_in, [[],[]] )
    ch_versions = ch_versions.mix(CONTIGS_FAIDX.out.versions)
    CONTIGS_FAIDX.out.fai.join( fasta_in )
                    .map{ meta, fai, fasta -> fasta }
                    .set{ scaf_ref }
    CONTIGS_FAIDX.out.fai.join( fasta_in )
                    .map{ meta, fai, fasta -> fai }
                    .set{ scaf_ref_fai }
    YAHS( bed_in, scaf_ref, scaf_ref_fai )
    ch_versions = ch_versions.mix(YAHS.out.versions)
    SCAFFOLDS_FAIDX(YAHS.out.scaffolds_fasta, [[],[]])
    ch_versions = ch_versions.mix(SCAFFOLDS_FAIDX.out.versions)
    bed_in.map{ meta, bed -> meta}.set{ch_meta}

    // Prepare contact pairs for cooler
    YAHS.out.binary.join(YAHS.out.scaffolds_agp)
                    .combine(scaf_ref)
                    .combine(scaf_ref_fai)
                    .map{meta, binary, agp, fa, fai -> [meta, binary, agp, fai]}
                    .set{ch_merge}
    JUICER_PRE(ch_merge)
    ch_versions = ch_versions.mix(JUICER_PRE.out.versions)
    
    // Bin contact pairs
    JUICER_PRE.out.pairs.join(bed_in)
                        .combine(Channel.of(cool_bin))
                        .set{ch_juicer}
    CHROM_SIZES(SCAFFOLDS_FAIDX.out.fai)
    COOLER_CLOAD(ch_juicer, CHROM_SIZES.out.chrom_sizes)
    ch_versions = ch_versions.mix(COOLER_CLOAD.out.versions)
    
    // Generate a multi-resolution cooler file by coarsening
    COOLER_CLOAD.out.cool.map{ meta, cools, cool_bin-> [meta, cools]}
                         .set{ch_cool}
    COOLER_ZOOMIFY(ch_cool)
    ch_versions = ch_versions.mix(COOLER_ZOOMIFY.out.versions)

    // Create contact map in pretext format
    SCAFFOLDS_FAIDX.out.fai.map{ meta, fai -> fai }.set{fai}
    PREPARE_PRETEXTMAP_INPUT(JUICER_PRE.out.pairs, fai)
    ch_versions = ch_versions.mix(PREPARE_PRETEXTMAP_INPUT.out.versions)
    PRETEXTMAP(PREPARE_PRETEXTMAP_INPUT.out.pairs, [])
    ch_versions = ch_versions.mix(PRETEXTMAP.out.versions)

    PRETEXTSNAPSHOT(PRETEXTMAP.out.pretext)
    ch_versions = ch_versions.mix(PRETEXTSNAPSHOT.out.versions)

    // Generate HiC Map
    JUICER_TOOLS_PRE(JUICER_PRE.out.pairs, CHROM_SIZES.out.chrom_sizes, 'yahs_scaffolds')
    ch_versions = ch_versions.mix(JUICER_TOOLS_PRE.out.versions)

    emit:
    alignments_sorted = JUICER_PRE.out.pairs
    fasta = YAHS.out.scaffolds_fasta
    chrom_sizes = CHROM_SIZES.out.chrom_sizes
    cool = COOLER_CLOAD.out.cool
    mcool = COOLER_ZOOMIFY.out.mcool
    snapshots = PRETEXTSNAPSHOT.out.image
    hic = JUICER_TOOLS_PRE.out.hic
    versions = ch_versions.ifEmpty(null)
}
