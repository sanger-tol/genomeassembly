process PRETEXT_MAP {
    tag "$meta.id"
    label 'process_medium'

    input:
    tuple val(meta), path(pairs)
    val(fai)

    output:
    tuple val(meta), path("*pretext"), emit: pretext
    path "versions.yml"           , emit: versions
    
    when:
    task.ext.when == null || task.ext.when


    script:
    def args = task.ext.args ?: ''
    """
    prepare_pretext.sh $fai $pairs > $meta.id".prep";
    cat $meta.id".prep" | PretextMap $args -o $meta.id".pretext";

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        PretextMap: \$(PretextMap | grep 'Version' | cut -f3 -d' ')
    END_VERSIONS
    """


}
