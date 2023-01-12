process JUICER_TOOLS_PRE {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::java-jdk=8.0.112" : null)
    container "${ workflow.containerEngine == 'singularity' &&
                    !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/java-jdk:8.0.112--1' :
        'quay.io/biocontainers/java-jdk:8.0.112--1' }"

    input:
    tuple val(meta), path(pairs)
    path sizes
    val prefix

    output:
    tuple val(meta), path("*hic"), emit: hic
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when


    script:
    def juicer_tools_jar = task.ext.juicer_tools_jar ?: ''
    def juicer_jvm_params = task.ext.juicer_jvm_params ?: ''
    """
    java ${juicer_jvm_params} -jar ${projectDir}/bin/${juicer_tools_jar} pre \\
        ${pairs} ${prefix}.hic ${sizes}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        juicer tools: \$(java ${juicer_jvm_params} -jar ${projectDir}/bin/${juicer_tools_jar} -V | grep "Juicer Tools Version" | sed 's/Juicer Tools Version //')
    END_VERSIONS
    """
}
