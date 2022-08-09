include { COOLER_CLOAD    } from '../../modules/nf-core/modules/cooler/cload/main.nf' 
include { COOLER_ZOOMIFY         } from '../../modules/nf-core/modules/cooler/zoomify/main.nf' 
include { SAMTOOLS_FAIDX } from '../../modules/nf-core/modules/samtools/faidx/main.nf'
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
    bed_in.map{ meta, bed -> meta}.set{ch_meta}

    // Prepare contact pairs for cooler
    YAHS.out.binary.join(YAHS.out.scaffolds_agp)
                    .combine(fasta_in)
                    .map{meta, binary, agp, fa, fai -> [meta, binary, agp, fai]}
                    .set{ch_merge}
    JUICER_PRE(ch_merge)
    ch_versions = ch_versions.mix(JUICER_PRE.out.versions)
    
    // Bin contact paits
    JUICER_PRE.out.pairs.join(bed_in)
                        .set{ch_juicer}
    CHROM_SIZES(SAMTOOLS_FAIDX.out.fai)
    COOLER_CLOAD(ch_juicer, cool_bin, CHROM_SIZES.out.chrom_sizes)
    ch_versions = ch_versions.mix(COOLER_CLOAD.out.versions)
    
    // Generate a multi-resolution cooler file by coarsening
    COOLER_CLOAD.out.cool.map{ meta, cool_bin, cools -> [meta, cools]}
                         .set{ch_cool}
    COOLER_ZOOMIFY(ch_cool)
    ch_versions = ch_versions.mix(COOLER_ZOOMIFY.out.versions)

    // Create contact map in pretext format
    SAMTOOLS_FAIDX.out.fai.map{ meta, fai -> fai }.set{fai}
    PRETEXT_MAP(JUICER_PRE.out.pairs, fai)
    ch_versions = ch_versions.mix(PRETEXT_MAP.out.versions)

    PRETEXT_SNAPSHOT(PRETEXT_MAP.out.pretext)
    ch_versions = ch_versions.mix(PRETEXT_SNAPSHOT.out.versions)

    // Generate HiC Map
    JUICER_TOOLS_PRE(JUICER_PRE.out.pairs, CHROM_SIZES.out.chrom_sizes, 'yahs_scaffolds')
    ch_versions = ch_versions.mix(JUICER_TOOLS_PRE.out.versions)

    JUICER_SNAPSHOT(JUICER_TOOLS_PRE.out.hic)
    ch_versions = ch_versions.mix(JUICER_SNAPSHOT.out.versions)

    emit:
    alignments_sorted = JUICER_PRE.out.pairs
    fasta = YAHS.out.scaffolds_fasta
    chrom_sizes = CHROM_SIZES.out.chrom_sizes
    stats = GFASTATS.out.stats
    cool = COOLER_CLOAD.out.cool
    mcool = COOLER_ZOOMIFY.out.mcool
    snapshots = PRETEXT_SNAPSHOT.out.snapshot
    hic = JUICER_TOOLS_PRE.out.hic
    png = JUICER_SNAPSHOT.out.png
    versions = ch_versions.ifEmpty(null)

}
