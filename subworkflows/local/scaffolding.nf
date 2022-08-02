//TODO: install as nf-core module
include { COOLER_CLOAD    } from '../../../modules/modules/cooler/cload/main.nf' 
include { COOLER_ZOOMIFY         } from '../../../modules/modules/cooler/zoomify/main.nf' 
//TODO: install as nf-core module  and synchronize with polishing
include { SAMTOOLS_FAIDX } from '../../../modules/modules/samtools/faidx/main.nf'
include { YAHS           } from '../../modules/local/yahs.nf'
include { JUICER_PRE         } from '../../modules/local/juicer_pre.nf'
include { JUICER_SNAPSHOT         } from '../../modules/local/juicer_snapshot.nf'
include { JUICER_TOOLS_PRE         } from '../../modules/local/juicer_tools_pre.nf'
include { PRETEXT_MAP         } from '../../modules/local/pretext_map.nf'
include { PRETEXT_SNAPSHOT         } from '../../modules/local/pretext_snapshot.nf'
include { CHROM_SIZES    } from '../../modules/local/chrom_sizes.nf'
include { GFASTATS       } from '../../modules/local/gfastats.nf'

workflow SCAFFOLDING {
    take:  
    bed_in  // tuple(meta, bed)
    fasta_in   // tuple (fasta, fai)
    ifbreak // val: if/break/contigs/or/not/for/yahs
    motif // val: restriction/enzyme/motif
    resolutions // val: resolution/parameter/for/yahs
    cool_bin  // val: cooler cload parameter 
    
    main:
    ch_versions = Channel.empty()
    YAHS( bed_in, fasta_in, ifbreak, motif, resolutions )
    ch_versions = ch_versions.mix(YAHS.out.versions)
    SAMTOOLS_FAIDX(YAHS.out.scaffolds_fasta)
    ch_versions = ch_versions.mix(SAMTOOLS_FAIDX.out.versions)
    GFASTATS(YAHS.out.scaffolds_fasta)
    ch_versions = ch_versions.mix(GFASTATS.out.versions)
    YAHS.out.binary.map{ meta, binary -> [ 'yahs', binary ] }
        .join(YAHS.out.scaffolds_agp.map {meta, agp -> [ 'yahs', agp]})
        .join(  fasta_in.map{ fa, fai -> [ 'yahs', fai]} )
        .join(bed_in.map{ meta, bed -> ['yahs', meta] })
        .map{ yahs, binary, agp, fai, meta -> [meta, binary, agp, fai]}
        .set{ch_merge}
    JUICER_PRE(ch_merge)
    ch_versions = ch_versions.mix(JUICER_PRE.out.versions)
    JUICER_PRE.out.pairs.map{meta, it -> ['juicer', it]}
                        .join( bed_in.map{ meta, bed -> ['juicer', meta] })
                        .map{ juicer, it, meta -> [meta, it, []]}
                        .set{ch_juicer}
    CHROM_SIZES(SAMTOOLS_FAIDX.out.fai)
    COOLER_CLOAD(ch_juicer, cool_bin, CHROM_SIZES.out.chrom_sizes)
    ch_versions = ch_versions.mix(COOLER_CLOAD.out.versions)
    COOLER_CLOAD.out.cool.map{ meta, cool_bin, cools -> ['cload', cools]}
                         .join( bed_in.map{ meta, bed -> ['cload', meta] })
                         .map{ cload, cools, meta -> [meta, cools]}
                         .set{ch_cool}
    COOLER_ZOOMIFY(ch_cool)
    ch_versions = ch_versions.mix(COOLER_ZOOMIFY.out.versions)
    SAMTOOLS_FAIDX.out.fai.map{ meta, fai -> fai }.set{fai}
    PRETEXT_MAP(JUICER_PRE.out.pairs, fai)
    ch_versions = ch_versions.mix(PRETEXT_MAP.out.versions)
    PRETEXT_SNAPSHOT(PRETEXT_MAP.out.pretext)
    ch_versions = ch_versions.mix(PRETEXT_SNAPSHOT.out.versions)
    JUICER_TOOLS_PRE(JUICER_PRE.out.pairs, CHROM_SIZES.out.chrom_sizes, 'yahs_scaffolds')
    ch_versions = ch_versions.mix(JUICER_TOOLS_PRE.out.versions)
    JUICER_SNAPSHOT(JUICER_TOOLS_PRE.out.hic)
    ch_versions = ch_versions.mix(JUICER_SNAPSHOT.out.versions)

    emit:
    fasta = YAHS.out.scaffolds_fasta
    chrom_sizes = CHROM_SIZES.out.chrom_sizes
    stats = GFASTATS.out.stats
    alignments_sorted = JUICER_PRE.out.pairs
    cool = COOLER_CLOAD.out.cool
    mcool = COOLER_ZOOMIFY.out.mcool
    snapshots = PRETEXT_SNAPSHOT.out.snapshot
    hic = JUICER_TOOLS_PRE.out.hic
    png = JUICER_SNAPSHOT.out.png
    versions = ch_versions.ifEmpty(null)

}
