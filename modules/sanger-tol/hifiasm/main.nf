process HIFIASM {
    tag "${meta.id}"
    label 'process_high_memory'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/37/373e515588f1ee1cc0eb4c68fe83db690ccd372d134036464bcb24a2ef060b27/data' :
        'community.wave.seqera.io/library/hifiasm_htslib_samtools_gawk:6a66c53cc0723e53' }"

    input:
    tuple val(meta) , path(long_reads)        , path(ul_reads)
    tuple val(meta2), path(hic_cram)
    tuple val(meta3), path(paternal_kmer_dump), path(maternal_kmer_dump)
    tuple val(meta4), path(bin_files)

    output:
    tuple val(meta), path("*.bin")        , emit: bin_files        , optional: true
    tuple val(meta), path("*.gfa.gz")     , emit: assembly_graphs  , optional: true
    tuple val(meta), path("*.fa")         , emit: assembly_fasta   , optional: true
    tuple val(meta), path("*.ec.fa.gz")   , emit: corrected_reads  , optional: true
    tuple val(meta), path("*.ovlp.paf.gz"), emit: read_overlaps    , optional: true
    tuple val(meta), path("*.log")        , emit: log
    path  "versions.yml"                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"

    def long_reads_sorted = long_reads instanceof List ? long_reads.sort{ file -> file.name } : long_reads
    def hic_cram_sorted = hic_cram instanceof List ? hic_cram.sort{ file -> file.name } : hic_cram
    def ul_reads_sorted = ul_reads instanceof List ? ul_reads.sort{ file -> file.name } : ul_reads
    def ultralong = ul_reads ? "--ul ${ul_reads_sorted}" : ""

    if([paternal_kmer_dump, maternal_kmer_dump].any() && hic_cram) {
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
    if(hic_cram) {
        input_hic1 = "--h1 <(for f in ${hic_cram_sorted}; do samtools cat \$f | samtools fastq -n -f0x40 -F0xB00; done)"
        input_hic2 = "--h2 <(for f in ${hic_cram_sorted}; do samtools cat \$f | samtools fastq -n -f0x80 -F0xB00; done)"
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

    for graph in *.gfa; do
        bn=\$(basename \${graph} .gfa)
        gawk '
            BEGIN { OFS = "\\t" }
            /^S/ {
                print ">" \$2
                print \$3
            }
        ' \$graph > \${bn}.fa
    done

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
    echo "" | bgzip > ${prefix}.r_utg.gfa.gz
    touch ${prefix}.r_utg.fa
    echo "" | bgzip > ${prefix}.bp.p_ctg.gfa.gz
    touch ${prefix}.bp.p_ctg.fa
    echo "" | bgzip > ${prefix}.p_utg.gfa.gz
    touch ${prefix}.p_utg.fa
    echo "" | bgzip > ${prefix}.p_ctg.gfa.gz
    touch ${prefix}.p_ctg.fa
    echo "" | bgzip > ${prefix}.a_ctg.gfa.gz
    touch ${prefix}.a_ctg.fa
    echo "" | bgzip > ${prefix}.hap1.p_ctg.gfa.gz
    touch ${prefix}.hap1.p_ctg.fa
    echo "" | bgzip > ${prefix}.hap2.p_ctg.gfa.gz
    touch ${prefix}.hap2.p_ctg.fa
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
