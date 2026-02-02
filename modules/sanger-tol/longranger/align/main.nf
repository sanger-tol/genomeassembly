process LONGRANGER_ALIGN {
    tag "${meta.id}"
    label 'process_high'

    // NOTE: As longranger is proprietary software, no container or conda environment
    // can be provided for this module. Please configure your own container and then
    // pass the path to this container as part of a Nextflow configuration file.

    input:
    tuple val(meta) , path(fastq, stageAs: "10X_inputs/*")
    tuple val(meta2), path(reference)

    output:
    tuple val(meta), path("${prefix}/outs/possorted_bam.bam")    , emit: bam
    tuple val(meta), path("${prefix}/outs/possorted_bam.bam.bai"), emit: bai
    tuple val(meta), path("${prefix}/outs/summary.csv")          , emit: csv
    path("versions.yml")                                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def local_setup = !(task.ext.args =~ /--jobmode=lsf|--jobmode=sge/) ? "--localcores=${task.cpus}" : ""
    def args        = task.ext.args   ?: ""
    prefix          = task.ext.prefix ?: "${meta.id}"
    """
    longranger align \\
        --id=${prefix} \\
        --fastqs=10X_inputs \\
        --reference=${reference} \\
        ${local_setup} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        longranger: \$(longranger align --version | grep longranger | sed 's/.*(//' | sed 's/).*//')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}/outs
    touch ${prefix}/outs/possorted_bam.bam
    touch ${prefix}/outs/possorted_bam.bam.bai
    touch ${prefix}/outs/summary.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        longranger: \$(longranger align --version | grep longranger | sed 's/.*(//' | sed 's/).*//')
    END_VERSIONS
    """
}
