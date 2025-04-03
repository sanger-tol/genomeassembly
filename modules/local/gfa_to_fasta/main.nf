process GFA_TO_FASTA {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/5c/5cb24e1a833f9f340973d1250ef6d7ca8ff32995aa4f1c1850e66ff0a991c613/data' :
        'community.wave.seqera.io/library/samtools_gawk:2f6ad5ac0a3fef78' }"

    input:
    tuple val(meta), path(gfa)

    output:
    tuple val(meta), path("*.fa.gz"), emit: fasta
    path "versions.yml"             , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    zcat ${gfa} |\
        awk '/^S/ {
            print ">"\$2
            print \$3
        }' |\
        bgzip -@ ${task.cpus} > ${prefix}.fa.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gawk: \$(echo \$(awk --version | head -n1 | cut -f3 -d' ' | sed 's/,//'))
        samtools: \$(samtools --version |& sed '1!d ; s/samtools //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}.fa.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gawk: \$(echo \$(awk --version | head -n1 | cut -f3 -d' ' | sed 's/,//'))
        samtools: \$(samtools --version |& sed '1!d ; s/samtools //')
    END_VERSIONS
    """
}
