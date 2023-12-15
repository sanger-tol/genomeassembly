process FASTA_CONCAT {
    label 'process_single'

    conda "conda-forge::gawk=5.1.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gawk:5.1.0' :
        'quay.io/biocontainers/gawk:5.1.0' }"

    input:
    tuple val(meta), path(primary), path(haplotigs)

    output:
    tuple val(meta), path("*fa"), emit: fasta_merged
    path "versions.yml"         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    cat $primary $haplotigs > primary_and_haplotigs.fa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        GNU coreutils: \$(cat --version | head -n 1 | awk '{print \$NF}')
    END_VERSIONS
    """
}
