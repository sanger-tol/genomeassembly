process GFASTATS {
    tag "${meta.id}"
    label 'process_single'

    conda (params.enable_conda ? "bioconda::gfastats=1.3.5" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gfastats:1.3.5--hd03093a_0' :
        'quay.io/biocontainers/gfastats:1.3.5--hd03093a_0' }"
    input:
    tuple val(meta), path(fasta)

    output:
    path("*stats"), emit: stats
    path "versions.yml",     emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    gfastats ${fasta} ${args} > ${fasta}.stats
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gfatstas: \$(gfastats --version | grep gfastats | cut -d' ' -f2)
    END_VERSIONS
    """
}
