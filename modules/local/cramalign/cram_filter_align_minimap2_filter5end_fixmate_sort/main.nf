process CRAM_FILTER_MINIMAP2_FILTER5END_FIXMATE_SORT {
    tag "$meta.id"
    label "process_high"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/65/65076cea537714273d691bebd9ad2be1a64d9f6fe5c992528389c8d4c6820d32/data' :
        'community.wave.seqera.io/library/minimap2_samtools_staden_io_lib_perl:50147563892d7138' }"

    input:
    tuple val(meta), path(cram), path(crai), val(chunkn), val(range), path(reference)

    output:
    tuple val(meta), path("*.bam"), emit: bam
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args1 = task.ext.args1 ?: ''
    def args2 = task.ext.args2 ?: ''
    def args3 = task.ext.args3 ?: ''
    def args4 = task.ext.args4 ?: ''
    def args5 = task.ext.args5 ?: ''
    def args6 = task.ext.args6 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    cram_filter ${args1} -n ${range[0]}-${range[1]} ${cram} - |\\
        samtools fastq ${args2} - |\\
        minimap2 -t${task.cpus} ${args3} ${reference} - |\\
        filter_five_end.pl |\\
        samtools fixmate ${args4} - - |\\
        samtools view -h ${args5} |\\
        samtools sort ${args6} -@${task.cpus} -T ${prefix}_tmp -o ${prefix}.bam -

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        staden_io_lib: \$(io_lib-config --version 2>&1)
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//' )
        minimap2: \$(minimap2 --version | sed 's/minimap2 //g')
    END_VERSIONS
    """

    stub:
    def prefix  = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        staden_io_lib: \$(io_lib-config --version 2>&1)
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//' )
        minimap2: \$(minimap2 --version | sed 's/minimap2 //g')
    END_VERSIONS
    """
}
