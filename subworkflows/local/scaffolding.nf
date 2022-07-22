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
    bed   // path: path/to/merged.bed
    fasta // path: path/to/contigs
    fasta_fai   // path: path/to/fasta/index
    ifbreak // val: if/break/contigs/or/not/for/yahs
    motif // val: restriction/enzyme/motif
    resolutions // val: resolution/parameter/for/yahs

    main:
    ch_versions = Channel.empty()
    YAHS( [ [id: 'yahs'], bed, fasta, fasta_fai, ifbreak, 
            motif, resolutions ] )
    SAMTOOLS_FAIDX(YAHS.out.scaffolds_fasta)
    GFASTATS(YAHS.out.scaffolds_fasta)
    CHROM_SIZES(SAMTOOLS_FAIDX.out.fai.collect{it[1]})
    ch_fasta_fai = Channel.of(['yahs',fasta_fai ])
    YAHS.out.binary.map{ meta, binary -> [ meta.id.toString(), binary] }
        .join(YAHS.out.scaffolds_agp.map {meta, agp -> [meta.id.toString(), agp]})
        .join( ch_fasta_fai )
        .map{ meta, binary, agp, fai -> [[id:'juicer'], binary, agp, fai] }
        .set{ch_merge}
    JUICER_PRE(ch_merge)
    JUICER_PRE.out.pairs.map{ it -> [[id: 'cooler'], it, []]}
                        .set{ch_juicer}
    COOLER_CLOAD(ch_juicer, 1000, CHROM_SIZES.out.chrom_sizes)
    COOLER_CLOAD.out.cool.map{ meta, cool_bin, cools -> [meta, cools]}
                         .set{ch_cool}
    COOLER_ZOOMIFY(ch_cool)
    //COOLER_CLOAD([[id:'yahs'], JUICER_PRE.out.pairs, fasta_fai], Channel.of(1000), CHROM_SIZES.out.chrom_sizes)
    emit:
    fasta = YAHS.out.scaffolds_fasta
}
