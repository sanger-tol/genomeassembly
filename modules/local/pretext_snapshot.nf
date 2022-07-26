process PRETEXT_SNAPSHOT {
    tag "$meta.id"
    label 'process_medium'

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
