process MITOHIFI_MITOHIFI {
    tag "$meta.id"
    label 'process_high'

    // Docker image available at the biocontainers Dockerhub
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://biocontainers/mitohifi:2.2_cv1':
        'docker.io/biocontainers/mitohifi:2.2_cv1' }"

    input:
    tuple val(meta), path(reads), path(contigs)
    path ref_fa
    path ref_gb
    val mito_code

    output:
    tuple val(meta), path("*fasta"), emit: fasta
    tuple val(meta), path("*gb"), emit: gb, optional: true
    tuple val(meta), path("*gff"), emit: gff, optional: true
    path  "versions.yml"                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    if (reads) {
    """
    mitohifi.py -r ${reads} -f ${ref_fa} -g ${ref_gb} -o ${mito_code} -t $task.cpus ${args}
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mitohifi: \$( mitohifi.py --version 2>&1 | head -n1 | sed 's/^.*MitoHiFi //; s/ .*\$//' )
    END_VERSIONS
    """
    }
    else if (contigs) {
    """
    mitohifi.py -c ${contigs} -f ${ref_fa} -g ${ref_gb} -o ${mito_code} -t $task.cpus ${args}
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mitohifi: \$( mitohifi.py --version 2>&1 | head -n1 | sed 's/^.*MitoHiFi //; s/ .*\$//' )
    END_VERSIONS
        """
    }
}
