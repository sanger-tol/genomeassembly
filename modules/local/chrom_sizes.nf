process CHROM_SIZES {
    label 'process_nompi'

    conda (params.enable_conda ? "conda-forge::gawk=5.1.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gawk:5.1.0' :
        'quay.io/biocontainers/gawk:5.1.0' }"

    input:
    path fai 
    
    output:
    path "*chrom.sizes", emit: chrom_sizes
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    chrom_sizes.sh $fai > $fai".chrom.sizes"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        GNU Awk: \$(echo \$(awk --version 2>&1) | grep -i awk | sed 's/GNU Awk //; s/,.*//')
    END_VERSIONS
    """
}
