process GET_CALCUTS_PARAMS {
    label 'process_single'

    conda (params.enable_conda ? "conda-forge::python" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3' :
        'python:3' }"

    input:
    tuple val(meta), path(model_fk)

    output:
    env cutoffs        , emit: cutoffs
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = meta.prefix ?: ''
    """
    cutoffs=`get_calcuts_params_from_model_fk.py $model_fk`

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        \$(python --version)
    END_VERSIONS
    """
}
