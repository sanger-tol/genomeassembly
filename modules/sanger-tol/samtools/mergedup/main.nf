process SAMTOOLS_MERGEDUP {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/8c/8c5d2818c8b9f58e1fba77ce219fdaf32087ae53e857c4a496402978af26e78c/data'
        : 'community.wave.seqera.io/library/htslib_samtools:1.23.1--5b6bb4ede7e612e5'}"

    input:
    tuple val(meta) , path(input_files, stageAs: "?/*"), path(index_files, stageAs: "?/*")
    tuple val(meta2), path(fasta), path(fai), path(gzi)

    output:
    tuple val(meta), path("*.bam")      , emit: bam,  optional: true
    tuple val(meta), path("*.cram")     , emit: cram, optional: true
    tuple val(meta), path("*.{csi,crai}"), emit: index, optional: true
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
        --threads ${task.cpus-1} \\
        ${args} \\
        - \\
        ${reference} \\
        ${input_files} |\\
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
