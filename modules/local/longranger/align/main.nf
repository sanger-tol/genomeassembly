process LONGRANGER_ALIGN {
    // To use in cluster mode, some extra configurations is needed.
    // Visit tests/modules/longranger/align/nextflow.config for an example.

    tag "$meta.id"
    label 'process_high'

    input:
    tuple val(meta) , path(reference)
    tuple val(meta2), path(fastqs, stageAs: "10X_inputs/*")

    output:
    tuple val(meta)    , path("${meta.id}/outs/possorted_bam.bam")    , emit: bam
    tuple val(meta)    , path("${meta.id}/outs/possorted_bam.bam.bai"), emit: bai
    tuple val(meta)    , path("${meta.id}/outs/summary.csv")          , emit: csv
    path("versions.yml")                                              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def local_setup = !(task.ext.args =~ /--jobmode=lsf|--jobmode=sge/) ? "--localcores=${task.cpus}" : ""
    def args        = task.ext.args   ?: ""
    def prefix      = task.ext.prefix ?: "${meta.id}"
    """
    longranger align \\
        --id=${prefix} \\
        --fastqs=10X_inputs \\
        --sample=${prefix} \\
        --reference=${reference} \\
        ${local_setup} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        longranger: \$(longranger align --version | grep longranger | sed 's/.*(//' | sed 's/).*//')
    END_VERSIONS
    """

    stub:
    def prefix      = reference - ~/\.fa/
    """
    mkdir -p refdata-${prefix}/outs
    touch refdata-${prefix}/outs/possorted_bam.bam
    touch refdata-${prefix}/outs/possorted_bam.bam.bai
    touch refdata-${prefix}/outs/summary.csv
    """
}
