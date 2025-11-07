process FASTXALIGN_PYFASTXINDEX {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/63/637ecb387eaafe0b1689e3e32c5eda589e016cfd46c482946425181f69f0733e/data' :
        'community.wave.seqera.io/library/htslib_minimap2_pyfastx_samtools_click:bfd8f60cc27aa6d6' }"

    input:
    tuple val(meta), path(fastx)

    output:
    tuple val(meta), path("*.fxi"), stdout, emit: index
    path "versions.yml"                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // WARNING: This module includes the slice_fasta.py script as a module binary in
    // ${moduleDir}/resources/usr/bin/slice_fasta.py. To use this module, you will
    // either have to copy this file to ${projectDir}/bin or set the option
    // nextflow.enable.moduleBinaries = true
    // in your nextflow.config file.
    def args       = task.ext.args  ?: ''
    """
    slice_fasta.py index \\
        ${fastx} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        slice_fasta.py: \$(slice_fasta.py --version)
    END_VERSIONS
    """

    stub:
    """
    touch ${fastx}.fxi
    ## output dummy count to stdout
    echo -n 100

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        slice_fasta.py: \$(slice_fasta.py --version)
    END_VERSIONS
    """
}
