process YAHS_MAKEPAIRSFILE {
    tag "${meta.id}"
    label "process_medium"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/4c/4c5d6ee67b551df7a0641fc25727b66f356f2d0ae0bf9fe0d730e567210e6d2f/data' :
        'community.wave.seqera.io/library/htslib_samtools_yahs_coreutils_gawk:c423bb8e262dc97d' }"

    input:
    tuple val(meta), path(scaffolds_fai), path(scaffolds_agp), path(contigs_fai), path(contigs_contacts)

    output:
    tuple val(meta), path("*.pairs.gz"), emit: pairs
    tuple val("${task.process}"), val('juicer_pre'), eval('juicer pre --version'), emit: versions_juicerpre, topic: versions
    tuple val("${task.process}"), val('gawk'), eval('awk -Wversion | sed "1!d; s/.*Awk //; s/,.*//"'), emit: versions_gawk, topic: versions
    tuple val("${task.process}"), val('sort'), eval('sort --version | grep -o -E "[0-9]+(\\.[0-9]+)+"'), emit: versions_sort, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def buffer_mem  = (task.memory.toGiga() / 2).round()
    """
    gawk '
        BEGIN {
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
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}.pairs.gz
    """
}
