process BEDTOOLS_BAMTOBEDSORT {
    tag "$meta.id"
    label "process_high"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/3e/3ed03d6b4fc208960b4053d5d0a320fc6ca1ef39eef09a89aab9f43790638db6/data' :
        'community.wave.seqera.io/library/bedtools_htslib_samtools_coreutils:a52bfee84c9f20c4' }"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.bed"), emit: sorted_bed
    tuple val("${task.process}"), val('bedtools'), eval('bedtools --version | sed -e "s/bedtools v//g"'), emit: versions_bedtools, topic: versions
    tuple val("${task.process}"), val('samtools'), eval('samtools version | sed "1!d;s/.* //"'), emit: versions_samtools, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix      = task.ext.prefix ?: "${meta.id}"
    def args        = task.ext.args   ?: "-hu"
    def args2       = task.ext.args2  ?: ""
    def args3       = task.ext.args3  ?: ""
    def st_cores    = task.cpus > 4 ? 4 : task.cpus
    def buffer_mem  = (task.memory.toGiga() / 2).round()
    """
    samtools view \\
        -@${st_cores} \\
        ${args} \\
        ${bam} | \\
    bamToBed ${args2} -i stdin | \\
    sort ${args3} \\
        --parallel=${task.cpus} \\
        -S ${buffer_mem}G \\
        -T . > \\
    ${prefix}.bed
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bed
    """
}
