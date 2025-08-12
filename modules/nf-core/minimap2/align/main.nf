process MINIMAP2_ALIGN {
    tag "$meta.id"
    label 'process_high'

    // Note: the versions here need to match the versions used in the mulled container below and minimap2/index
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/66/66dc96eff11ab80dfd5c044e9b3425f52d818847b9c074794cf0c02bfa781661/data' :
        'community.wave.seqera.io/library/minimap2_samtools:33bb43c18d22e29c' }"

    input:
    tuple val(meta), path(reads)
    tuple val(meta2), path(reference)
    val bam_format
    val bam_index_extension
    val cigar_paf_format
    val cigar_bam
    val output_meta

    output:
    tuple val(meta_out), path("*.paf")                       , optional: true, emit: paf
    tuple val(meta_out), path("*.bam")                       , optional: true, emit: bam
    tuple val(meta_out), path("*.bam.${bam_index_extension}"), optional: true, emit: index
    path "versions.yml"                                      , emit: versions
    when:
    task.ext.when == null || task.ext.when

    script:
    def args  = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def args3 = task.ext.args3 ?: ''
    def args4 = task.ext.args4 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def bam_index = bam_index_extension ? "${prefix}.bam##idx##${prefix}.bam.${bam_index_extension} --write-index" : "${prefix}.bam"
    def bam_output = bam_format ? "-a | samtools sort -@ ${task.cpus-1} -o ${bam_index} ${args2}" : "-o ${prefix}.paf"
    def cigar_paf = cigar_paf_format && !bam_format ? "-c" : ''
    def set_cigar_bam = cigar_bam && bam_format ? "-L" : ''
    def bam_input = "${reads.extension}".matches('sam|bam|cram')
    def samtools_reset_fastq = bam_input ? "samtools reset --threads ${task.cpus-1} $args3 $reads | samtools fastq --threads ${task.cpus-1} $args4 |" : ''
    def query = bam_input ? "-" : reads
    def target = reference ?: (bam_input ? error("BAM input requires reference") : reads)
    if (!(output_meta in ['reads', 'reference', 'combine'])){
        log.error("Error: Invalid value for output_meta")
    }
    if(output_meta == 'reads')     { meta_out = meta         }
    if(output_meta == 'reference') { meta_out = meta2        }
    if(output_meta == 'combine'  ) { meta_out = meta + meta2 }
    """
    $samtools_reset_fastq \\
    minimap2 \\
        $args \\
        -t $task.cpus \\
        $target \\
        $query \\
        $cigar_paf \\
        $set_cigar_bam \\
        $bam_output


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2: \$(minimap2 --version 2>&1)
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def output_file = bam_format ? "${prefix}.bam" : "${prefix}.paf"
    def bam_index = bam_index_extension ? "touch ${prefix}.bam.${bam_index_extension}" : ""
    def bam_input = "${reads.extension}".matches('sam|bam|cram')
    def target = reference ?: (bam_input ? error("BAM input requires reference") : reads)

    """
    touch $output_file
    ${bam_index}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2: \$(minimap2 --version 2>&1)
    END_VERSIONS
    """
}
