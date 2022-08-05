#!/usr/bin/env nextflow

nextflow.enable.dsl=2
include { SCAFFOLDING } from './subworkflows/local/scaffolding'


workflow {
    def meta = [:]
    meta.id = "idPhyRufi1"
//    bed="/lustre/scratch123/tol/projects/darwin/data/insects/Physocephala_rufipes/working/idPhyRufi1.hifiasm.20220613/scaff/merge.mkdup.bed"
    bed = "/lustre/scratch123/tol/projects/darwin/data/insects/Physocephala_rufipes/working/idPhyRufi1.hifiasm.20220613/scaffolding/merge.mkdup.bed"  
    ch_bed = Channel.of([meta, bed])

//    fasta = "/lustre/scratch123/tol/projects/darwin/data/insects/Physocephala_rufipes/working/idPhyRufi1.hifiasm.20220613/polishing/purged.consensus.fa"
    fasta = "/lustre/scratch123/tol/projects/darwin/data/insects/Physocephala_rufipes/working/idPhyRufi1.hifiasm.20220613/scaffolding/ref/ref.fa"

    ch_fa = Channel.of([fasta, fasta+".fai"])

    SCAFFOLDING (ch_bed, ch_fa, true, '', '', 1000)
}
