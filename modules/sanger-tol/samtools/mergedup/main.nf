process SAMTOOLS_MERGEDUP {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/e5/e5598451c6d348cce36191bafe1911ad71e440137d7a329da946f2b0dbb0e7f3/data' :
        'community.wave.seqera.io/library/htslib_samtools:1.23--cde2c40a51d6f752' }"

    input:
    tuple val(meta) , path(input)
    tuple val(meta2), path(fasta), path(fai),  path(gzi)

    output:
    tuple val(meta), path("*.bam")      , emit: bam,  optional: true
    tuple val(meta), path("*.cram")     , emit: cram, optional: true
    tuple val(meta), path("*.csi")      , emit: csi,  optional: true
    tuple val(meta), path("*.crai")     , emit: crai, optional: true
    tuple val(meta), path("*.metrics")  , emit: metrics
    tuple val("${task.process}"), val('samtools'), eval('samtools version | sed "1!d;s/.* //"'), emit: versions_samtools, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args      = task.ext.args  ?: ''
    def args2     = task.ext.args2 ?: ''
    def prefix    = task.ext.prefix ?: "${meta.id}"
    def reference = fasta ? "--reference ${fasta}" : ""
    def extension = args2.contains("--output-fmt sam") ? "sam" :
                    args2.contains("--output-fmt cram") ? "cram" :
                    "bam"
    """
    samtools merge \\
        ${args} \\
        - \\
        ${input} |\\
    samtools markdup \\
        -T ${prefix} \\
        -f ${prefix}.metrics \\
        --threads ${task.cpus} \\
        $reference \\
        $args2 \\
        - \\
        ${prefix}.${extension}
    """

    stub:
    def args2     = task.ext.args2 ?: ''
    def prefix    = task.ext.prefix ?: "${meta.id}"
    def extension = args2.contains("--output-fmt sam") ? "sam" :
                    args2.contains("--output-fmt cram") ? "cram" :
                    "bam"
    """
    touch ${prefix}.${extension}
    touch ${prefix}.metrics
    """
}
