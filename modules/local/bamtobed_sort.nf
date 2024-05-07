//
// Copied from https://github.com/sanger-tol/treeval/blob/28309b7a1faf3aee5627f497c7cfa62d12ac65b8/modules/local/bamtobed_sort.nf
// from Sanger TOL treeval pipeline
//


process BAMTOBED_SORT {
    tag "$meta.id"
    label "process_high"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-9d3a458f6420e5712103ae2af82c94d26d63f059:60b54b43045e8cf39ba307fd683c69d4c57240ce-0' :
        'biocontainers/mulled-v2-9d3a458f6420e5712103ae2af82c94d26d63f059:60b54b43045e8cf39ba307fd683c69d4c57240ce-0' }"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.bed"), emit: sorted_bed
    path "versions.yml"           , emit: versions

    script:
    def prefix      = args.ext.prefix ?: "${meta.id}"
    def st_cores    = task.cpus > 4 ? 4 : "${task.cpus}"
    def buffer_mem  = (task.memory.toGiga() / 2).round()
    """
    samtools view -@${st_cores} -u -F0x400 ${bam} | bamToBed | sort -k4 --parallel=${task.cpus} -S ${buffer_mem}G -T . > ${prefix}_merged_sorted.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
        bedtools: \$(bedtools --version | sed -e "s/bedtools v//g")
    END_VERSIONS
    """

    stub:
    def prefix = args.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_merged_sorted.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
        bedtools: \$(bedtools --version | sed -e "s/bedtools v//g")
    END_VERSIONS
    """
}
