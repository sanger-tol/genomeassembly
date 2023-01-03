process PRETEXT_SNAPSHOT {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::pretextsnapshot:0.0.4" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pretextsnapshot:0.0.4--h9f5acd7_1' :
        'quay.io/biocontainers/pretextsnapshot:0.0.4--h9f5acd7_1' }"

    input:
    tuple val(meta), path(pretext)

    output:
    tuple val(meta), path("*snapshots/*png"),  emit: snapshot
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when


    script:
    def args = task.ext.args ?: ''
    """
    PretextSnapshot -m $pretext $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        PretextSnapshot: \$(PretextSnapshot | grep 'Version' | cut -f3 -d' ')
    END_VERSIONS
    """


}
