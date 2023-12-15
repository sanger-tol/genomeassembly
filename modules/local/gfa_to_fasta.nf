process GFA_TO_FASTA {
    tag "$meta.id"
    label 'process_high'

    conda "conda-forge::gawk=5.1.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gawk:5.1.0' :
        'quay.io/biocontainers/gawk:5.1.0' }"

    input:
    tuple val(meta), path(gfa)

    output:
    tuple val(meta), path("*.fa"),         emit: fasta
    path "versions.yml",                   emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    awk '/^S/{print ">"\$2;print \$3}' $gfa > ${prefix}.fa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        GNU Awk: \$(echo \$(awk --version | head -n1 | cut -f3 -d' ' | sed 's/,//'))
    END_VERSIONS
    """
}
