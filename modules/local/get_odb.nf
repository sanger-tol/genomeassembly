// -----
// Version from Sanger Genomenote pipeline by @priyanka-surana
// https://github.com/sanger-tol/genomenote/blob/383f23e6b7a89f9aad6b85c8f7320b5c5825de73/modules/local/get_odb.nf
// ------


process GET_ODB {
    tag "${meta.id}"
    label 'process_single'

    conda (params.enable_conda ? "conda-forge::requests=2.26.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/requests:2.26.0' :
        'quay.io/biocontainers/requests:2.26.0' }"
    input:
    tuple val(meta), val(tolid)

    output:
    path("*.busco_odb.csv"), emit: csv
    path "versions.yml",     emit: versions

    script:
    """
    get_odb.py ${tolid} ${tolid}.busco_odb.csv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        get_odb.py: \$(get_odb.py --version | cut -d' ' -f2)
    END_VERSIONS
    """
}
