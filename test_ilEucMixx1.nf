#!/usr/bin/env nextflow

nextflow.enable.dsl=2
include { SCAFFOLDING } from './subworkflows/local/scaffolding'


workflow {
    def meta = [:]
    meta.id = "ilEucMixx1"
//    bed="/lustre/scratch123/tol/projects/darwin/data/insects/Euclidia_mi/working/ilEucMixx1.hifiasm.20220324/scaffolding/merge.mkdup.bed"
    bed="/lustre/scratch123/tol/projects/darwin/data/insects/Euclidia_mi/working/ilEucMixx1.hifiasm.20220324/scaffolding/merge.mkdup.bed"
    
    ch_bed = Channel.of([meta, bed])

//    fasta = "/lustre/scratch123/tol/projects/darwin/data/insects/Euclidia_mi/working/ilEucMixx1.hifiasm.20220324/polishing/purged.consensus.fa"
    fasta = "/lustre/scratch123/tol/projects/darwin/data/insects/Euclidia_mi/working/ilEucMixx1.hifiasm.20220324/scaffolding/ref/ref.fa"
    ch_fa = Channel.of([fasta, fasta+".fai"])

    
    SCAFFOLDING (ch_bed, ch_fa, true, '', '', 1000)
}
