process HIFIASM {
    tag "${meta.id}"
    label 'process_high_memory'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/17/17053ea7b1ae932805a745499662dfbca1c2267597a0126ef48d33bc90bc1d6a/data'
        : 'community.wave.seqera.io/library/hifiasm_htslib_samtools_findutils_gawk:fb65b90c39bfe422'}"

    input:
    tuple val(meta), path(long_reads, arity: '1..*'), path(ul_reads, arity: '0..*')
    tuple val(meta2), path(hic_cram, arity: '0..*')
    tuple val(meta3), path(paternal_kmer_dump), path(maternal_kmer_dump)
    tuple val(meta4), path(bin_files)
    tuple val(meta5), path(previous_log, stageAs: "log/*")

    output:
    tuple val(meta), path("*.bin"), emit: bin_files, optional: true
    tuple val(meta), path("*.gfa.gz"), emit: assembly_graphs, optional: true
    tuple val(meta), path("*.fa"), emit: assembly_fasta, optional: true
    tuple val(meta), path("*.ec.fa.gz"), emit: corrected_reads, optional: true
    tuple val(meta), path("*.ovlp.paf.gz"), emit: read_overlaps, optional: true
    tuple val(meta), path("*.bed.gz"), emit: bed, optional: true
    tuple val(meta), path("*.log"), emit: log
    tuple val("${task.process}"), val('hifiasm'), eval("hifiasm --version"), emit: versions_hifiasm, topic: versions
    tuple val("${task.process}"), val('samtools'), eval('samtools --version | head -1 | sed -e "s/samtools //"'), emit: versions_samtools, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"

    // Sort all input files for reproducibility
    def long_reads_sorted = long_reads.sort { file -> file.name }
    def hic_cram_sorted = hic_cram.sort { file -> file.name }
    def ul_reads_sorted = ul_reads.sort { file -> file.name }

    // Check if we are trying to do both trio assembly and phased assembly
    if ([paternal_kmer_dump, maternal_kmer_dump].any() && hic_cram) {
        log.error("ERROR: hifiasm trio binning mode and Hi-C phasing can not be used at the same time.")
    }

    // Validate trio inputs and set up input trio arg
    def input_trio = ""
    if ([paternal_kmer_dump, maternal_kmer_dump].any()) {
        if (![paternal_kmer_dump, maternal_kmer_dump].every()) {
            log.error("ERROR: Either the maternal or paternal kmer dump is missing!")
        }
        else {
            input_trio = "-1 ${paternal_kmer_dump} -2 ${maternal_kmer_dump}"
        }
    }

    // Set up hic input arg
    def input_hic1 = ""
    def input_hic2 = ""
    if (hic_cram) {
        input_hic1 = "--h1 <(samtools cat ${hic_cram_sorted} | samtools fastq -n -f0x40 -F0xB00)"
        input_hic2 = "--h2 <(samtools cat ${hic_cram_sorted} | samtools fastq -n -f0x80 -F0xB00)"
    }

    // Set up ultralong reads input
    def ultralong = ul_reads ? "--ul <(cat ${ul_reads_sorted})" : ""

    // Configure log input so that logs are always written to a new file
    def copy_previous_log = previous_log ? "cat log/* > ${prefix}.log" : ""
    """
    ${copy_previous_log}

    hifiasm \\
        ${args} \\
        -t ${task.cpus} \\
        ${input_trio} \\
        ${input_hic1} \\
        ${input_hic2} \\
        ${ultralong} \\
        -o ${prefix} \\
        ${long_reads_sorted} \\
        2> >( tee -a ${prefix}.log >&2 )

    if [ -f ${prefix}.ec.fa ]; then
        bgzip -@${task.cpus} ${prefix}.ec.fa
    fi

    if [ -f ${prefix}.ovlp.paf ]; then
        bgzip -@${task.cpus} ${prefix}.ovlp.paf
    fi

    ## Note there is presently a bug in Hifiasm when running with --n-hap >=3,
    ## where all contigs get the "h1tg" schema instead of the actual haplotype
    ## number. In the name of simpler code we just arrange a fix for all Hi-C
    ## GFA files.
    ##
    ## https://github.com/chhylp123/hifiasm/issues/489
    find . -maxdepth 1 -name "*.gfa" | while read graph; do
        if [[ "\$graph" =~ hic\\.hap([0-9]+)\\.p_ctg\\.gfa ]]; then
            hap_num="\${BASH_REMATCH[1]}"
            bn=\$(basename "\$graph" .gfa)

            gawk -v i="\$hap_num" '{ gsub("h1tg", "h" i "tg"); print }' "\$graph" > "\${bn}.fixed.gfa"
            mv \${bn}.fixed.gfa \${bn}.gfa
        fi
    done

    ## Convert GFA to FASTA.
    find . -maxdepth 1 -name "*.gfa" | while read graph; do
        bn=\$(basename "\$graph" .gfa)
        gawk 'BEGIN { OFS = "\\t" } /^S/ { print ">" \$2; print \$3 }' "\$graph" > "\${bn}.fa"
    done

    ## gzip all GFA and BED output files
    find . -name "*.gfa" -exec bgzip -@${task.cpus} {} \\;
    find . -name "*.bed" -exec bgzip -@${task.cpus} {} \\;
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.ec.bin
    touch ${prefix}.ovlp.source.bin
    touch ${prefix}.ovlp.reverse.bin
    touch ${prefix}.hic.tlb.bin
    touch ${prefix}.hic.lk.bin
    for f in r_utg bp.p_ctg p_utg p_ctg a_ctg hap1.p_ctg hap2.p_ctg
    do
      touch ${prefix}.\${f}.gfa
      bgzip ${prefix}.\${f}.gfa
      touch ${prefix}.\${f}.fa
    done
    touch ${prefix}.ec.fa
    bgzip ${prefix}.ec.fa
    touch ${prefix}.ovlp.paf
    bgzip ${prefix}.ovlp.paf
    touch ${prefix}.log
    """
}
