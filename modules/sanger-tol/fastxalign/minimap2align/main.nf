process FASTXALIGN_MINIMAP2ALIGN {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/e3/e3e96ee4627b21402c6503d365bdc27e61f3ef5ae6dbeeda02b70a7157f2de9d/data' :
        'community.wave.seqera.io/library/htslib_minimap2_pyfastx_samtools_pruned:58888de833fa75ea' }"

    input:
    tuple val(meta),  path(fastx), path(fxi)
    tuple val(meta2), path(index), path(reference)
    tuple val(chunkn), val(range)
    val bam_format

    output:
    tuple val(meta), path("*.bam")   , emit: bam, optional: true
    tuple val(meta), path("*.paf.gz"), emit: paf, optional: true
    tuple val("${task.process}"), val('slice_fasta.py'), eval('slice_fasta.py --version'), emit: versions_slice_fasta, topic: versions
    tuple val("${task.process}"), val('samtools'), eval('samtools version | sed "1!d;s/.* //"'), emit: versions_samtools, topic: versions
    tuple val("${task.process}"), val('minimap2'), eval('minimap2 --version | sed "s/minimap2 //g"'), emit: versions_minimap2, topic: versions

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
        minimap2 -t${task.cpus} ${args1} ${index} - \\
        ${bam_output}
    """

    stub:
    def prefix  = task.ext.prefix ?: "${fastx}.${chunkn}.${meta.id}"
    """
    touch ${prefix}.bam
    echo "" | gzip > ${prefix}.paf.gz
    """
}
