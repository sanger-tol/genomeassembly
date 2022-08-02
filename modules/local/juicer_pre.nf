process JUICER_PRE {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::yahs:1.2a.1--h7132678_0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?        'https://depot.galaxyproject.org/singularity/yahs:1.2a.1--h7132678_0' :
        'quay.io/biocontainers/yahs:1.2a.1--h7132678_0' }"

    input:
    tuple val(meta), path(binary), path(agp), path(fai)

    output:
    tuple val(meta), path("*alignments_sorted.txt"), emit: pairs
    path "versions.yml"           , emit: versions
    
    when:
    task.ext.when == null || task.ext.when


    script:
    def args = task.ext.args ?: ''
    """
    juicer pre ${binary} ${agp} ${fai} ${args} | LC_ALL=C sort -k2,2d -k6,6d \\
        --parallel=8 -S50G | awk '\$3>=0 && \$7>=0'> alignments_sorted.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        juicer: \$(juicer 2>&1 | grep -w "Version:" | cut -f2 -d' ')
    END_VERSIONS
    """
}
