process CHROM_SIZES {
    tag "$meta.id"
    label 'process_nompi'

    conda "conda-forge::gawk=5.1.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gawk:5.1.0' :
        'quay.io/biocontainers/gawk:5.1.0' }"

    input:
    tuple val(meta), path(fai)

    output:
    path "*chrom.sizes", emit: chrom_sizes
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    cut -f1-2 $fai > \$(basename $fai).chrom.sizes

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gawk: \$(echo \$(awk --version 2>&1) | grep -i awk | sed 's/GNU Awk //; s/,.*//')
    END_VERSIONS
    """
}
