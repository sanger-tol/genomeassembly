/*
Copied from
https://github.com/NBISweden/Earth-Biogenome-Project-pilot/blob/5ec2002638055bb8396857a8ee418bf86188fc59/modules/local/purgedups/pbcstat.nf
*/
process PURGEDUPS_PBCSTAT {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::purge_dups=1.2.6" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/purge_dups:1.2.6--h7132678_0':
        'quay.io/biocontainers/purge_dups:1.2.6--h7132678_0' }"

    input:
    tuple val(meta), path(paf_alignment)

    output:
    tuple val(meta), path("*.PB.stat")    , emit: stat
    tuple val(meta), path("*.PB.base.cov"), emit: basecov
    path "versions.yml"                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    pbcstat \\
        $args \\
        $paf_alignment

    for PBFILE in PB.*; do mv \$PBFILE ${prefix}.\$PBFILE; done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        purgedups: \$( purge_dups -h |& sed '3!d; s/.*: //' )
    END_VERSIONS
    """
}
