process MAKE_PAIRS_FILE {
    tag "${meta.id}"
    label "process_low"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/ac/acf7567da73ff97edaab983f28f434826d242cc5684429724f23ad0b0ce99f5f/data' :
        'community.wave.seqera.io/library/samtools_yahs_gawk:0e36292a5b98b34c' }"

    input:
    tuple val(meta), path(scaffolds_fai), path(scaffolds_agp), path(contigs_fai), path(contigs_bin)

    output:
    tuple val(meta), path("*.pairs.gz"), emit: pairs
    path "versions.yml"                , emit: versions

    script:
    def prefix        = task.ext.prefix ?: "${meta.id}"
    """
    gawk '
        BEGIN {
            OFS = " "
            print "## pairs format v1.0"
        }
        { print "#chromsize:", \$1, \$2 }
        END {
            print "#columns:", "readID", "chr1", "pos1", "chr2", "pos2", "strand1", "strand2"
        }' ${scaffolds_fai} > ${prefix}.pairs

    juicer pre \\
        ${contigs_bin} \\
        ${scaffolds_agp} \\
        ${contigs_fai} |\\
        gawk '
            BEGIN { OFS = "\t" }
            \$3 > 0 && \$7 > 0 {print ".", \$2, \$3, \$6, \$7, ".", "."}
        ' |\\
        LC_ALL=C sort -k2,2d -k4,4d -S50G >>\\
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
