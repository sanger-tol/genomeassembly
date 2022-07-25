include { COOLER_CLOAD    } from '../../../modules/modules/cooler/cload/main.nf' 
include { COOLER_ZOOMIFY         } from '../../../modules/modules/cooler/zoomify/main.nf' 
//TODO: install as nf-core module  and synchronize with polishing
include { SAMTOOLS_FAIDX } from '../../../modules/modules/samtools/faidx/main.nf'
include { YAHS           } from '../../modules/local/yahs.nf'
include { JUICER_PRE         } from '../../modules/local/juicer_pre.nf'
include { CHROM_SIZES    } from '../../modules/local/chrom_sizes.nf'
include { GFASTATS       } from '../../modules/local/gfastats.nf'

workflow SCAFFOLDING {
    take:  
    bed_in  // tuple(meta, bed)
    fasta_in   // tuple (fasta, fai)
    ifbreak // val: if/break/contigs/or/not/for/yahs
    motif // val: restriction/enzyme/motif
    resolutions // val: resolution/parameter/for/yahs

    main:
    ch_versions = Channel.empty()
    YAHS( bed_in, fasta_in, ifbreak, motif, resolutions )
    SAMTOOLS_FAIDX(YAHS.out.scaffolds_fasta)
    GFASTATS(YAHS.out.scaffolds_fasta)
    YAHS.out.binary.map{ meta, binary -> [ 'yahs', binary ] }
        .join(YAHS.out.scaffolds_agp.map {meta, agp -> [ 'yahs', agp]})
        .join(  fasta_in.map{ fa, fai -> [ 'yahs', fai]} )
        .join(bed_in.map{ meta, bed -> ['yahs', meta] })
        .map{ yahs, binary, agp, fai, meta -> [meta, binary, agp, fai]}
        .set{ch_merge}
    JUICER_PRE(ch_merge)
    JUICER_PRE.out.pairs.map{meta, it -> ['juicer', it]}
                        .join( bed_in.map{ meta, bed -> ['juicer', meta] })
                        .map{ juicer, it, meta -> [meta, it, []]}
                        .set{ch_juicer}
    CHROM_SIZES(SAMTOOLS_FAIDX.out.fai)
    COOLER_CLOAD(ch_juicer, 1000, CHROM_SIZES.out.chrom_sizes)
    COOLER_CLOAD.out.cool.map{ meta, cool_bin, cools -> ['cload', cools]}
                         .join( bed_in.map{ meta, bed -> ['cload', meta] })
                         .map{ cload, cools, meta -> [meta, cools]}
                         .set{ch_cool}
    COOLER_ZOOMIFY(ch_cool)

    emit:
    fasta = YAHS.out.scaffolds_fasta
    chrom_sizes = CHROM_SIZES.out.chrom_sizes
    stats = GFASTATS.out.stats
    alignments_sorted = JUICER_PRE.out.pairs
    cool = COOLER_CLOAD.out.cool
    mcool = COOLER_ZOOMIFY.out.mcool
}
