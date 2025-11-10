process FASTXALIGN_MINIMAP2ALIGN {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/63/637ecb387eaafe0b1689e3e32c5eda589e016cfd46c482946425181f69f0733e/data' :
        'community.wave.seqera.io/library/htslib_minimap2_pyfastx_samtools_click:bfd8f60cc27aa6d6' }"

    input:
    tuple val(meta), path(fastx), path(fxi), val(chunkn), val(range), path(reference)
    val bam_format

    output:
    tuple val(meta), path("*.bam")   , emit: bam, optional: true
    tuple val(meta), path("*.paf.gz"), emit: paf, optional: true
    path "versions.yml"              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // WARNING: This module includes the slice_fasta.py script as a module binary in
    // ${moduleDir}/resources/usr/bin/slice_fasta.py. To use this module, you will
    // either have to copy this file to ${projectDir}/bin or set the option
    // nextflow.enable.moduleBinaries = true
    // in your nextflow.config file.
    def args1       = task.ext.args1  ?: ''
    def args2       = task.ext.args2  ?: ''
    def args3       = task.ext.args3  ?: ''
    def prefix      = task.ext.prefix ?: "${fastx}.${chunkn}.${meta.id}"
    def post_filter = args2 ? "samtools view -h ${args2} - |" : ''
    def sort_bam    = "samtools sort -@ ${task.cpus > 1 ? task.cpus - 1 : 1} -o ${prefix}.bam -T ${prefix}_sort_tmp ${args3} -"
    def bam_output  = bam_format      ? "-a | ${post_filter} ${sort_bam}" : "| bgzip -@ ${task.cpus} > ${prefix}.paf.gz"
    """
    slice_fasta.py slice ${fastx} ${range[0]} ${range[1]} | \\
        minimap2 -t${task.cpus} ${args1} ${reference} - \\
        ${bam_output}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        slice_fasta.py: \$(slice_fasta.py --version)
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//' )
        minimap2: \$(minimap2 --version | sed 's/minimap2 //g')
    END_VERSIONS
    """

    stub:
    def prefix  = task.ext.prefix ?: "${fastx}.${chunkn}.${meta.id}"
    """
    touch ${prefix}.bam
    echo "" | gzip > ${prefix}.paf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        slice_fasta.py: \$(slice_fasta.py --version)
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//' )
        minimap2: \$(minimap2 --version | sed 's/minimap2 //g')
    END_VERSIONS
    """
}
