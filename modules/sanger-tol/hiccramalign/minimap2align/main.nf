process HICCRAMALIGN_MINIMAP2ALIGN {
    tag "$meta.id"
    label "process_high"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/65/65858e733832166824cfd05291fc456bdf219b02baa3944c2c92efad86a6ee7f/data' :
        'community.wave.seqera.io/library/htslib_minimap2_samtools_gawk_perl:6729620c63652154' }"

    input:
    tuple val(meta), path(cram), path(crai), val(chunkn), val(range), path(reference)

    output:
    tuple val(meta), path("*.bam"), emit: bam
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // WARNING: This module includes the filter_five_end.pl script as a module binary in
    // ${moduleDir}/resources/usr/bin/filter_five_end.pl. To use this module, you will
    // either have to copy this file to ${projectDir}/bin or set the option
    // nextflow.enable.moduleBinaries = true
    // in your nextflow.config file.
    def args1 = task.ext.args1 ?: ''
    def args2 = task.ext.args2 ?: ''
    def args3 = task.ext.args3 ?: ''
    def args4 = task.ext.args4 ?: ''
    def args5 = task.ext.args5 ?: ''
    def args6 = task.ext.args6 ?: ''
    def prefix  = task.ext.prefix ?: "${cram}.${chunkn}.${meta.id}"
    """
    samtools cat ${args1} -r "#:${range[0]}-${range[1]}" ${cram} |\\
        samtools fastq ${args2} - |\\
        minimap2 -t${task.cpus} ${args3} ${reference} - |\\
        gawk '
            \$1 ~ /^\\@/ { print \$0 }
            \$1 !~ /^\\@/ && and(\$2, 64) > 0 { print 1 \$0 }
            \$1 !~ /^\\@/ && and(\$2, 64) == 0 { print 2 \$0 }
        ' |\\
        filter_five_end.pl |\\
        gawk '
            BEGIN { OFS="\\t" }
            \$1 ~ /^\\@/ { print \$0 }
            \$1 !~ /^\\@/ { \$2 = and(\$2, compl(2048)); print substr(\$0, 2) }
        ' |\\
        samtools fixmate ${args4} - - |\\
        samtools view -h ${args5} |\\
        samtools sort ${args6} -@${task.cpus} -T ${prefix}_tmp -o ${prefix}.bam -

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//' )
        minimap2: \$(minimap2 --version | sed 's/minimap2 //g')
        gawk: \$(gawk --version | grep -o -E "[0-9]+(\\.[0-9]+)+" | head -n1)
    END_VERSIONS
    """

    stub:
    def prefix  = task.ext.prefix ?: "${cram}.${chunkn}.${meta.id}"
    """
    touch ${prefix}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//' )
        minimap2: \$(minimap2 --version | sed 's/minimap2 //g')
        gawk: \$(gawk --version | grep -o -E "[0-9]+(\\.[0-9]+)+" | head -n1)
    END_VERSIONS
    """
}
