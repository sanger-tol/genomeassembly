process KEEP_SEQNAMES {
    label 'process_single'

    conda "conda-forge::coreutils=9.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
        'ubuntu:20.04' }"

    input:
    tuple val(meta), path(fa)

    output:
    path "*seq.lst"    , emit: seqlist
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = meta.prefix ?: ''
    """
    grep '>' $fa | cut -f1 | sed 's/>//' > seq.lst
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cut: \$(cut --version | head -n 1 | awk '{print \$NF}')
    END_VERSIONS
    """
}
