
process HAPMAKER {
    tag "$meta.id"
    label 'process_low'

    container 'ghcr.io/nbisweden/fastk_genescopefk_merquryfk:1.2'

    input:
    tuple val(meta), path(matktab)
    tuple val(meta2), path(patktab)
    tuple val(meta3), path(childktab)

    output:
    tuple val(meta), path("*mat*.hap.ktab*", hidden: true ), emit: mathap_ktab
    tuple val(meta), path("*pat*.hap.ktab*", hidden: true), emit: pathap_ktab
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def mat_ktab = matktab ? "${matktab.find{ it.toString().endsWith(".ktab") }}" : ''
    def pat_ktab = patktab ? "${patktab.find{ it.toString().endsWith(".ktab") }}" : ''
    def child_ktab = childktab ? "${childktab.find{ it.toString().endsWith(".ktab") }}" : ''
    def HAPmaker_VERSION = '1.2'
    """
    HAPmaker $args -T$task.cpus ${mat_ktab} ${pat_ktab} ${child_ktab} 2>hapmaker.log
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        HAPmaker: $HAPmaker_VERSION
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_mat.hap.ktab
    touch ${prefix}_pat.hap.ktab
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        HAPmaker: $HAPmaker_VERSION
    END_VERSIONS
    """
}
