#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { PREPARE_INPUT } from './subworkflows/local/prepare_input'

workflow {
    yaml = "/lustre/scratch124/tol/projects/darwin/users/kk16/development/nextflow/sanger-tol-genomeassembly/idPhyRufi1_test.yaml"

    ch_yaml = Channel.of(yaml)

    PREPARE_INPUT(ch_yaml)
}
