process BED_CHUNKS {
    label 'process_single'

    conda (params.enable_conda ? "conda-forge::gawk=5.1.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gawk:5.1.0' :
        'quay.io/biocontainers/gawk:5.1.0' }"

    input:
    tuple val(meta), path(fai)
    val chunks

    output:
    path "*bed"        , emit: coords
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    bed_chunks.sh $fai $chunks

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        GNU Awk: \$(echo \$(awk --version 2>&1) | grep -i awk | sed 's/GNU Awk //; s/,.*//')
        bed_chunks.sh: \$(bed_chunks.sh | tail -n 1 | cut -d' ' -f2)
    END_VERSIONS
    """
}
