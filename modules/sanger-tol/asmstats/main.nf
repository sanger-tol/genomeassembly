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
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // WARNING: This module includes the asmstats script as a module binary in
    // ${moduleDir}/resources/usr/bin/asmstats. To use this module, you will
    // either have to copy this file to ${projectDir}/bin or set the option
    // nextflow.enable.moduleBinaries = true in your nextflow.config file.
    def args    = task.ext.args ?: ''
    def prefix  = task.ext.prefix ?: "${assembly.getName()}"
    def VERSION = "1.0.0"
    """
    asmstats \\
        $args \\
        ${assembly} \\
        > ${prefix}.stats

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        asmstats: ${VERSION}
        seqtk: \$(seqtk |& sed '/Version/!d; s/.* //')
    END_VERSIONS
    """

    stub:
    def args    = task.ext.args ?: ''
    def prefix  = task.ext.prefix ?: "${meta.id}"
    def VERSION = "1.0.0"
    """
    touch ${prefix}.stats

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        asmstats: \$(asmstats --version)
        seqtk: \$(seqtk |& sed '/Version/!d; s/.* //')
    END_VERSIONS
    """
}
