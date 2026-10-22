process BEDCHUNKS_CREATE {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gawk:5.3.0' :
        'biocontainers/gawk:5.3.0' }"

    input:
    tuple val(meta), path(fai)
    val chunk_size

    output:
    tuple val(meta), path("*.bed"), emit: bed
    tuple val("${task.process}"), val('gawk'), eval('gawk -Wversion | sed "1!d; s/.*Awk //; s/,.*//"'), emit: versions_gawk, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix    = task.ext.prefix ?: "${meta.id}"
    """
    gawk -v prefix=${prefix} -v chunk_size=${chunk_size} \\
        'BEGIN { OFS = "\\t" }
        {
            chunk_id = NR % chunk_size
            print \$1, 0, \$2 >> prefix "." chunk_id ".bed"
        }' ${fai}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.1.bed
    touch ${prefix}.2.bed
    """
}
