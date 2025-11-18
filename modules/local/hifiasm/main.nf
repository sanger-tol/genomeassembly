process HIFIASM {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/70/709c9b732af5bc95b8f02346f374b8d9ef815087e30a969410d5e385625bf7a9/data' :
        'community.wave.seqera.io/library/hifiasm_samtools:7f8aef13225935f4' }"

    input:
    tuple val(meta) , path(long_reads)        , path(ul_reads)
    tuple val(meta1), path(paternal_kmer_dump), path(maternal_kmer_dump)
    tuple val(meta2), path(hic_reads)
    tuple val(meta3), path(bin_files)

    output:
    tuple val(meta), path("*.r_utg.gfa.gz")                             , emit: raw_unitigs      , optional: true
    tuple val(meta), path("*.bin")                                      , emit: bin_files        , optional: true
    tuple val(meta), path("*.p_utg.gfa.gz")                             , emit: processed_unitigs, optional: true
    tuple val(meta), path("${prefix}.{p_ctg,bp.p_ctg,hic.p_ctg}.gfa.gz"), emit: primary_contigs  , optional: true
    tuple val(meta), path("${prefix}.{a_ctg,hic.a_ctg}.gfa.gz")         , emit: alternate_contigs, optional: true
    tuple val(meta), path("${prefix}.*.hap1.p_ctg.gfa.gz")              , emit: hap1_contigs     , optional: true
    tuple val(meta), path("${prefix}.*.hap2.p_ctg.gfa.gz")              , emit: hap2_contigs     , optional: true
    tuple val(meta), path("*.ec.fa.gz")                                 , emit: corrected_reads  , optional: true
    tuple val(meta), path("*.ovlp.paf.gz")                              , emit: read_overlaps    , optional: true
    tuple val(meta), path("${prefix}.log")                              , emit: log
    path  "versions.yml"                                                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"

    def long_reads_sorted = long_reads instanceof List ? long_reads.sort{ file -> file.name } : long_reads
    def ul_reads_sorted = ul_reads instanceof List ? ul_reads.sort{ file -> file.name } : ul_reads
    def ultralong = ul_reads ? "--ul ${ul_reads_sorted}" : ""

    if([paternal_kmer_dump, maternal_kmer_dump].any() && hic_reads) {
        log.error("ERROR: hifiasm trio binning mode and Hi-C phasing can not be used at the same time.")
    }

    def input_trio = ""
    if([paternal_kmer_dump, maternal_kmer_dump].any()) {
        if(![paternal_kmer_dump, maternal_kmer_dump].every()) {
            log.error("ERROR: Either the maternal or paternal kmer dump is missing!")
        } else {
            input_trio = "-1 ${paternal_kmer_dump} -2 ${maternal_kmer_dump}"
        }
    }

    def input_hic1 = ""
    def input_hic2 = ""
    if(hic_reads) {
        input_hic1 = "--h1 <(for f in ${hic_reads}; do samtools cat \$f | samtools fastq -n -f0x40 -F0xB00; done)"
        input_hic2 = "--h2 <(for f in ${hic_reads}; do samtools cat \$f | samtools fastq -n -f0x80 -F0xB00; done)"
    }
    """
    hifiasm \\
        $args \\
        -t ${task.cpus} \\
        ${input_trio} \\
        ${input_hic1} \\
        ${input_hic2} \\
        ${ultralong} \\
        -o ${prefix} \\
        ${long_reads_sorted} \\
        2> >( tee ${prefix}.log >&2 )

    if [ -f ${prefix}.ec.fa ]; then
        bgzip -@${task.cpus} ${prefix}.ec.fa
    fi

    if [ -f ${prefix}.ovlp.paf ]; then
        bgzip -@${task.cpus} ${prefix}.ovlp.paf
    fi

    ## gzip all GFA output files
    find . -name "*.gfa" -exec bgzip -@${task.cpus} {} \\;

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hifiasm: \$(hifiasm --version 2>&1)
        samtools: \$(samtools --version |& sed '1!d ; s/samtools //')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.ec.bin
    touch ${prefix}.ovlp.source.bin
    touch ${prefix}.ovlp.reverse.bin
    touch ${prefix}.hic.tlb.bin
    touch ${prefix}.hic.lk.bin
    echo "" | bgzip ${prefix}.r_utg.gfa.gz
    echo "" | bgzip ${prefix}.bp.p_ctg.gfa.gz
    echo "" | bgzip ${prefix}.p_utg.gfa.gz
    echo "" | bgzip ${prefix}.p_ctg.gfa.gz
    echo "" | bgzip ${prefix}.a_ctg.gfa.gz
    echo "" | bgzip ${prefix}.hap1.p_ctg.gfa.gz
    echo "" | bgzip ${prefix}.hap2.p_ctg.gfa.gz
    echo "" | bgzip > ${prefix}.ec.fa.gz
    echo "" | bgzip > ${prefix}.ovlp.paf.gz
    touch ${prefix}.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hifiasm: \$(hifiasm --version 2>&1)
        samtools: \$(samtools --version |& sed '1!d ; s/samtools //')
    END_VERSIONS
    """
}
