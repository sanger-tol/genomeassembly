process ASMSTATS {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/c4/c47c4b4f55ce62f79edd9b4c5a26b7b1ed60684452dbbbd03f0be0d56c7b00be/data':
        'community.wave.seqera.io/library/seqtk_perl:37201934bb74266e' }"

    input:
    tuple val(meta), path(assembly)

    output:
    tuple val(meta), path("*.stats"), emit: stats
    // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    tuple val("${task.process}"), val('asmstats'), eval('echo 1.0.0'), emit: versions_asmstats, topic: versions
    tuple val("${task.process}"), val('seqtk'), eval('seqtk |& sed "/Version/!d; s/.* //"'), emit: versions_seqtk, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // WARNING: This module includes the asmstats script as a module binary in
    // ${moduleDir}/resources/usr/bin/asmstats. To use this module, you will
    // either have to copy this file to ${projectDir}/bin or set the option
    // nextflow.enable.moduleBinaries = true in your nextflow.config file.
    def args    = task.ext.args ?: ''
    def prefix  = task.ext.prefix ?: "${assembly.getName()}"
    """
    asmstats \\
        $args \\
        ${assembly} \\
        > ${prefix}.stats
    """

    stub:
    def prefix  = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.stats
    """
}
