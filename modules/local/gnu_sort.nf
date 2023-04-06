//
// Based on https://github.com/sanger-tol/genomenote/blob/2b17b1fb489178756e92ea4311edc1efd87d31c1/modules/local/gnu_sort.nf
// from Sanger genomenote pipeline by @priyanka-surana
//

process GNU_SORT {
    tag "$meta.id"
    label 'process_high'

    conda "conda-forge::sed=4.7" 
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
        'ubuntu:20.04' }"

    input:
    tuple val(meta), path(bed)

    output:
    tuple val(meta), path("*.sorted.bed"), emit: bed
    path "versions.yml",                   emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def mem = task.memory.giga
    """
    gnu_sort.sh "$args" ${task.cpus} ${mem}G $bed ${prefix}.sorted.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        GNU Sort: \$(echo \$(sort --version | grep sort | sed 's/.* //g'))
    END_VERSIONS
    """
}
