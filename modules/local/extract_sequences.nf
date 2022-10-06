process EXTRACT_SEQUENCES {
    label 'process_nompi'

    conda (params.enable_conda ? "conda-forge::seqtk=1.3" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/seqtk:1.3--h5bf99c6_3' :
        'quay.io/biocontainers/seqtk:1.3--h7132678_4' }"

    input:
    tuple val(meta), path(fasta) 
    path(fai)
    
    output:
    tuple val(meta), path("*subsequences.fasta"), emit: subseq
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    cat $fai | cut -f1 >  subseqs.lst
    seqtk subseq $fasta subseqs.lst > subsequences.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqtk: \$(seqtk 2>&1 | grep Version | cut -f2 -d' ')
    END_VERSIONS
    """
}
