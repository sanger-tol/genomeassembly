process PREPARE_PRETEXTMAP_INPUT {
    tag "$meta.id"
    label 'process_medium'

    conda "conda-forge::gawk=5.1.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gawk:5.1.0' :
        'quay.io/biocontainers/gawk:5.1.0' }"

    input:
    tuple val(meta), path(pairs)
    val(fai)

    output:
    tuple val(meta), path("*pairs.gz"), emit: pairs
    path "versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when


    script:
    def args = task.ext.args ?: ''
    """
    prepare_pretext.sh $fai $pairs | gzip -c > $meta.id".pairs.gz"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        GNU Awk: \$(echo \$(awk --version 2>&1) | grep -i awk | sed 's/GNU Awk //; s/,.*//')
    END_VERSIONS
    """


}
