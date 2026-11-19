process FASTXALIGN_PYFASTXINDEX {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/e3/e3e96ee4627b21402c6503d365bdc27e61f3ef5ae6dbeeda02b70a7157f2de9d/data' :
        'community.wave.seqera.io/library/htslib_minimap2_pyfastx_samtools_pruned:58888de833fa75ea' }"

    input:
    tuple val(meta), path(fastx)

    output:
    tuple val(meta), path(fastx), path("*.fxi"), stdout, emit: index
    tuple val("${task.process}"), val('slice_fasta.py'), eval('slice_fasta.py --version'), emit: versions_slice_fasta, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // WARNING: This module includes the slice_fasta.py script as a module binary in
    // ${moduleDir}/resources/usr/bin/slice_fasta.py. To use this module, you will
    // either have to copy this file to ${projectDir}/bin or set the option
    // nextflow.enable.moduleBinaries = true
    // in your nextflow.config file.
    def args       = task.ext.args  ?: ''
    """
    slice_fasta.py index \\
        ${fastx} \\
        ${args}
    """

    stub:
    """
    touch ${fastx}.fxi
    ## output dummy count to stdout
    echo -n 100
    """
}
