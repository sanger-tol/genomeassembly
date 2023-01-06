process JUICER_SNAPSHOT {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::java-jdk=8.0.112" : null)
    container "${ workflow.containerEngine == 'singularity' &&
                    !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/java-jdk:8.0.112--1' :
        'quay.io/biocontainers/java-jdk:8.0.112--1' }"


    input:
    tuple val(meta), path(hic)

    output:
    tuple val(meta), path("*png"), emit: png
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when


    script:
    def juicer_snapshot_jar = task.ext.juicer_snapshot_jar ?: ''
    def juicer_jvm_params = task.ext.juicer_jvm_params ?: ''
    """
    java ${juicer_jvm_params} -jar ${projectDir}/bin/${juicer_snapshot_jar} ${hic} \$(basename $hic)".png"  $args
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        JuicerSnapshot: \$(echo ${projectDir}/bin/${juicer_snapshot_jar} | sed 's/.*-//g' | sed 's/.jar//')
    END_VERSIONS
    """
}
