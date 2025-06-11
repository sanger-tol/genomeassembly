process BAMTOBED_SORT {
    tag "$meta.id"
    label "process_high"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/78/78feb7cf247c725d82ce8dc3b4cd2ea99be7e0034934485b4a44f79ee6edaf6a/data' :
        'community.wave.seqera.io/library/bedtools_samtools:73f3ef465b3bd116' }"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.bed"), emit: sorted_bed
    path "versions.yml"           , emit: versions

    script:
    def prefix      = task.ext.prefix ?: "${meta.id}"
    def args1       = task.ext.args1  ?: ""
    def args2       = task.ext.args2  ?: ""
    def st_cores    = task.cpus > 4 ? 4 : "${task.cpus}"
    def buffer_mem  = (task.memory.toGiga() / 2).round()
    """
    samtools view \\
        -@${st_cores} \\
        ${args1} ${bam} | \\
        bamToBed ${args2} -i stdin | \\
        sort -k4 --parallel=${task.cpus} -S ${buffer_mem}G -T . > \\
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
