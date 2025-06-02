process LONGRANGER_MKREF {
    tag "$meta.id"
    label 'process_medium'

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("refdata-*"), emit: folder
    path "versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    longranger mkref ${fasta}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        longranger: \$(longranger mkref --version | grep longranger | sed 's/.*(//' | sed 's/).*//')
    END_VERSIONS
    """

    stub:
    def prefix = fasta - ~/\.fa/
    """
    mkdir refdata-${prefix}
    mkdir -p refdata-${prefix}/fasta
    mkdir -p refdata-${prefix}/genes
    mkdir -p refdata-${prefix}/regions
    mkdir -p refdata-${prefix}/snps

    touch refdata-${prefix}/fasta/genome.fa
    touch refdata-${prefix}/fasta/genome.fa.amb
    touch refdata-${prefix}/fasta/genome.fa.ann
    touch refdata-${prefix}/fasta/genome.fa.bwt
    touch refdata-${prefix}/fasta/genome.fa.fai
    touch refdata-${prefix}/fasta/genome.fa.flat
    touch refdata-${prefix}/fasta/genome.fa.gdx
    touch refdata-${prefix}/fasta/genome.fa.pac
    touch refdata-${prefix}/fasta/genome.fa.sa

    touch refdata-${prefix}/genome

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        longranger: \$(longranger mkref --version | grep longranger | sed 's/.*(//' | sed 's/).*//')
    END_VERSIONS
    """
}
