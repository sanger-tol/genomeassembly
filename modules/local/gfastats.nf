process GFASTATS {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::gfastats=1.3.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gfastats:1.3.1--hd03093a_0' :
        'quay.io/biocontainers/gfastats:1.3.1--hd03093a_0' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*stats"), emit: stats
    path "versions.yml"           , emit: versions
    
    when:
    task.ext.when == null || task.ext.when


    script:
    def args = task.ext.args ?: ''
    """
    gfastats ${args} ${fasta} > ${fasta}.stats

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gfastats: \$(gfastats -v | head -n1 | cut -f2 -d' ')
    END_VERSIONS
    """


}
