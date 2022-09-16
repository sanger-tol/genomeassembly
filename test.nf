#!/usr/bin/env nextflow

include { POLISHING } from './subworkflows/local/polishing'

nextflow.enable.dsl=2

workflow {
    def meta = [:]
    meta.id = "idPhyRufi1"
    bam = "/lustre/scratch124/tol/projects/darwin/users/kk16/development/nextflow/freebayes-consensus-wf-dev/test.bam"
   ch_bam = Channel.of([meta, bam, bam+".bai"])

    fasta = "/lustre/scratch124/tol/projects/darwin/users/kk16/development/nextflow/freebayes-consensus-wf-dev/test.fa"

    ch_fa = Channel.of([fasta, fasta+".fai"])

    summary = "/lustre/scratch123/tol/projects/darwin/data/insects/Physocephala_rufipes/working/idPhyRufi1.hifiasm.20220613/polishing/idPhyRufi1/outs/summary.csv"
    
//    ch_summary = Channel.of(summary)


    groups = 2
//    ch_groups = Channel.of(groups)
    POLISHING (ch_bam, ch_fa, groups, summary)
//    POLISHING (bam, bam+".bai", fasta, fasta+".fai", groups, '')
}
