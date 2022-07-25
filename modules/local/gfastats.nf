process GFASTATS {
    tag "$meta.id"
    label 'process_medium'

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*stats"), emit: stats
    path "versions.yml"           , emit: versions
    
    when:
    task.ext.when == null || task.ext.when


    script:
    def args = task.ext.args ?: ''
    """
    gfastats ${args} ${fasta} > ${fasta}.stats

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gfastats: \$(gfastats -v | head -n1 | cut -f2 -d' ')
    END_VERSIONS
    """


}
