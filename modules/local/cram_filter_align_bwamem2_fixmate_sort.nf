//
// Copied from https://github.com/sanger-tol/treeval/blob/80554a803903183613d49690d5770eeadb3c42c9/modules/local/cram_filter_align_bwamem2_fixmate_sort.nf
// from Sanger TOL treeval pipeline
//

process CRAM_FILTER_ALIGN_BWAMEM2_FIXMATE_SORT {
    tag "$meta.id"
    label "process_high"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-50d89b457e04ed90fa0cbf8ebc3ae1b9ffbc836b:caf993da1689e8d42f5e4c113ffc9ef81d26df96-0' :
        'quay.io/biocontainers/mulled-v2-50d89b457e04ed90fa0cbf8ebc3ae1b9ffbc836b:caf993da1689e8d42f5e4c113ffc9ef81d26df96-0' }"

    input:
    tuple val(meta), path(cramfile), path(cramindex), val(from), val(to), val(base), val(chunkid), val(rglines), val(bwaprefix)

    output:
    tuple val(meta), path("*.bam"), emit: mappedbam
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    cram_filter -n ${from}-${to} ${cramfile} - | \\
        samtools fastq -F0xB00 -nt - | \\
        bwa-mem2 mem -p ${bwaprefix} -t${task.cpus} -5SPCp -H'${rglines}' - | \\
        samtools fixmate -mpu - - | \\
        samtools sort --write-index -l1 -@${task.cpus} -T ${base}_${chunkid}_sort_tmp -o ${prefix}_${base}_${chunkid}_mem.bam -

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//' )
        bwa-mem2: \$(bwa-mem2 --version | sed 's/bwa-mem2 //g')
    END_VERSIONS
    """
    // temp removal staden_io_lib: \$(echo \$(staden_io_lib --version 2>&1) | sed 's/^.*staden_io_lib //; s/Using.*\$//')

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_${base}_${chunkid}_mem.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//' )
        bwa-mem2: \$(bwa-mem2 --version | sed 's/bwa-mem2 //g')
    END_VERSIONS
    """
}
