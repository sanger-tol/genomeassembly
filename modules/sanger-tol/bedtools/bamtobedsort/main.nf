process BEDTOOLS_BAMTOBEDSORT {
    tag "$meta.id"
    label "process_high"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/cf/cf65e22485bc417775bef295e324b214d5a4ddfea5c3cbfedf8623bf8af55612/data' :
        'community.wave.seqera.io/library/bedtools_samtools_coreutils:43e34cbcf7cfff84' }"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.bed"), emit: sorted_bed
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix      = task.ext.prefix ?: "${meta.id}"
    def args1       = task.ext.args1  ?: ""
    def args2       = task.ext.args2  ?: ""
    def args3       = task.ext.args3  ?: ""
    def st_cores    = task.cpus > 4 ? 4 : "${task.cpus}"
    def buffer_mem  = (task.memory.toGiga() / 2).round()
    """
    samtools view \\
        -@${st_cores} \\
        ${args1} \\
        ${bam} | \\
    bamToBed ${args2} -i stdin | \\
    sort ${args3} \\
        --parallel=${task.cpus} \\
        -S ${buffer_mem}G \\
        -T . > \\
    ${prefix}.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
        bedtools: \$(bedtools --version | sed -e "s/bedtools v//g")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
        bedtools: \$(bedtools --version | sed -e "s/bedtools v//g")
    END_VERSIONS
    """
}
