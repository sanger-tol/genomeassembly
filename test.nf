#!/usr/bin/env nextflow

nextflow.enable.dsl=2
include { SCAFFOLDING } from './subworkflows/local/scaffolding'


workflow {
    def meta = [:]
//    meta.id = "idPhyRufi1"
    bed = "/lustre/scratch124/tol/projects/darwin/users/kk16/development/nextflow/yahs-wf/test_data/merge.mkdup.bed"
    ch_bed = Channel.of(bed)
    
 //   ch_bam = Channel.of([meta, bed])

    fasta = "/lustre/scratch124/tol/projects/darwin/users/kk16/development/nextflow/yahs-wf/test_data/test.fa"

    ch_fa = Channel.of(fasta)

    summary = "/lustre/scratch123/tol/projects/darwin/data/insects/Physocephala_rufipes/working/idPhyRufi1.hifiasm.20220613/polishing/idPhyRufi1/outs/summary.csv"
    
    ch_summary = Channel.of(summary)

    ch_fai = Channel.of(fasta+".fai")

    SCAFFOLDING (bed, fasta, fasta+".fai", true, '', '')
}
