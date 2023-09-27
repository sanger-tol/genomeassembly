//
// Copied from https://github.com/sanger-tol/treeval/blob/28309b7a1faf3aee5627f497c7cfa62d12ac65b8/modules/local/cram_filter_align_bwamem2_fixmate_sort.nf
// from Sanger TOL treeval pipeline
//

process CRAM_FILTER_ALIGN_BWAMEM2_FIXMATE_SORT {
    tag "$meta.id"
    label "process_high"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-50d89b457e04ed90fa0cbf8ebc3ae1b9ffbc836b:caf993da1689e8d42f5e4c113ffc9ef81d26df96-0' :
        'biocontainers/mulled-v2-50d89b457e04ed90fa0cbf8ebc3ae1b9ffbc836b:caf993da1689e8d42f5e4c113ffc9ef81d26df96-0' }"

    input:
    tuple val(meta), path(cramfile), path(cramindex), val(from), val(to), val(base), val(chunkid), val(rglines), val(bwaprefix)

    output:
    tuple val(meta), path("*.bam"), emit: mappedbam
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args1 = task.ext.args1 ?: ''
    def args2 = task.ext.args2 ?: ''
    def args3 = task.ext.args3 ?: ''
    def args4 = task.ext.args4 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // Please be aware one of the tools here required mem = 28 * reference size!!!
    """
    cram_filter -n ${from}-${to} ${cramfile} - | \\
        samtools fastq ${args1} | \\
        bwa-mem2 mem -p ${bwaprefix} -t${task.cpus} -5SPCp -H'${rglines}' - | \\
        samtools fixmate ${args3} - - | \\
        samtools sort ${args4} -@${task.cpus} -T ${base}_${chunkid}_sort_tmp -o ${prefix}_${base}_${chunkid}_mem.bam -

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//' )
        bwa-mem2: \$(bwa-mem2 --version | sed 's/bwa-mem2 //g')
    END_VERSIONS
    """
    // temp removal staden_io_lib: \$(echo \$(staden_io_lib --version 2>&1) | sed 's/^.*staden_io_lib //; s/Using.*\$//') CAUSES ERROR

    stub:
    def prefix  = task.ext.prefix ?: "${meta.id}"
    def base    = "45022_3#2"
    def chunkid = "1"
    """
    touch ${prefix}_${base}_${chunkid}_mem.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//' )
        bwamem2: \$(echo \$(bwa-mem2 version 2>&1) | sed 's/.* //')
    END_VERSIONS
    """
}
