process YAHS_MAKEPAIRSFILE {
    tag "${meta.id}"
    label "process_medium"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/dd/ddcb4312a2df1f479502951d398f1bb983f82728446784f8fa9c4e91074ffdb7/data' :
        'community.wave.seqera.io/library/samtools_yahs_coreutils_gawk:4f936f9704189766' }"

    input:
    tuple val(meta), path(scaffolds_fai), path(scaffolds_agp), path(contigs_fai), path(contigs_contacts)

    output:
    tuple val(meta), path("*.pairs.gz"), emit: pairs
    path "versions.yml"                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def buffer_mem  = (task.memory.toGiga() / 2).round()
    """
    gawk '
        BEGIN {
            OFS = "\t"
            print "## pairs format v1.0"
        }
        { print "#chromsize:", \$1, \$2 }
        END {
            print "#columns:", "readID", "chr1", "pos1", "chr2", "pos2", "strand1", "strand2"
        }' ${scaffolds_fai} > ${prefix}.pairs

    juicer pre \\
        ${contigs_contacts} \\
        ${scaffolds_agp} \\
        ${contigs_fai} |\\
    gawk '
        BEGIN { OFS = "\t" }
        \$3 > 0 && \$7 > 0 {print ".", \$2, \$3, \$6, \$7, ".", "."}
    ' |\\
    LC_ALL=C sort -k2,2d -k4,4d \\
        -S${buffer_mem}G \\
        -T . \\
        --parallel=${task.cpus} >>\\
        ${prefix}.pairs

    bgzip -@${task.cpus} ${prefix}.pairs

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        juicer_pre: \$(juicer_pre --version)
        gawk: \$(awk -Wversion | sed '1!d; s/.*Awk //; s/,.*//')
        sort: \$(sort --version | sed '1!d' | grep -o -E "[0-9]+(\\.[0-9]+)+")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}.pairs.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        juicer_pre: \$(juicer_pre --version)
        gawk: \$(awk -Wversion | sed '1!d; s/.*Awk //; s/,.*//')
        sort: \$(sort --version | sed '1!d' | grep -o -E "[0-9]+(\\.[0-9]+)+")
    END_VERSIONS
    """
}
