process HIFIASM {
    tag "$meta.id"
    label 'process_high'

    container "quay.io/sanger-tol/hifiasm_samtools:0.01"

    input:
    tuple val(meta), path(reads)
    path  paternal_kmer_dump
    path  maternal_kmer_dump
    path  hic_read1
    path  hic_read2
    path  hic_reads_cram

    output:
    tuple val(meta), path("*.r_utg.gfa")       , emit: raw_unitigs
    tuple val(meta), path("*.ec.bin")          , emit: corrected_reads
    tuple val(meta), path("*.ovlp.source.bin") , emit: source_overlaps
    tuple val(meta), path("*.ovlp.reverse.bin"), emit: reverse_overlaps
    tuple val(meta), path("*.bp.p_ctg.gfa")    , emit: processed_contigs, optional: true
    tuple val(meta), path("*.p_utg.gfa")       , emit: processed_unitigs, optional: true
    tuple val(meta), path("*.asm.p_ctg.gfa")   , emit: primary_contigs  , optional: true
    tuple val(meta), path("*.asm.a_ctg.gfa")   , emit: alternate_contigs, optional: true
    tuple val(meta), path("*.asm.hic.hap1.p_ctg.gfa")   , emit: hap1_contigs  , optional: true
    tuple val(meta), path("*.asm.hic.hap2.p_ctg.gfa")   , emit: hap2_contigs  , optional: true
    tuple val(meta), path("*.asm.hic.p_ctg.gfa")   , emit: hic_primary_contigs  , optional: true
    tuple val(meta), path("*.asm.hic.a_ctg.gfa")   , emit: hic_alternate_contigs  , optional: true
    tuple val(meta), path("*.asm.dip.hap1.p_ctg.gfa")  , emit: paternal_contigs , optional: true
    tuple val(meta), path("*.asm.dip.hap2.p_ctg.gfa")  , emit: maternal_contigs , optional: true
    tuple val(meta), path("*.log")             , emit: log
    path  "versions.yml"                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        exit 1, "This version of HIFIASM module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def hic_read1 = hic_reads_cram ? "for f in $hic_reads_cram; do samtools cat \$f | samtools fastq -n -f0x40 -F0xB00; done" : ""
    def hic_read2 = hic_reads_cram ? "for f in $hic_reads_cram; do samtools cat \$f | samtools fastq -n -f0x80 -F0xB00; done" : ""
    if ((paternal_kmer_dump) && (maternal_kmer_dump) && (hic_read1) && (hic_read2)) {
        error "Hifiasm Trio-binning and Hi-C integrated should not be used at the same time"
    } else if ((paternal_kmer_dump) && !(maternal_kmer_dump)) {
        error "Hifiasm Trio-binning requires maternal data"
    } else if (!(paternal_kmer_dump) && (maternal_kmer_dump)) {
        error "Hifiasm Trio-binning requires paternal data"
    } else if ((paternal_kmer_dump) && (maternal_kmer_dump)) {
        """
        hifiasm \\
            $args \\
            -o ${prefix}.asm \\
            -t $task.cpus \\
            -1 $paternal_kmer_dump \\
            -2 $maternal_kmer_dump \\
            $reads \\
            2> >( tee ${prefix}.stderr.log >&2 )


        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            hifiasm: \$(hifiasm --version 2>&1)
        END_VERSIONS
        """
    } else if ((hic_read1) && !(hic_read2)) {
        error "Hifiasm Hi-C integrated requires paired-end data (only R1 specified here)"
    } else if (!(hic_read1) && (hic_read2)) {
        error "Hifiasm Hi-C integrated requires paired-end data (only R2 specified here)"
    } else if ((hic_read1) && (hic_read2)) {
        """
        hifiasm \\
            $args \\
            -o ${prefix}.asm \\
            -t $task.cpus \\
            --h1 <($hic_read1) \\
            --h2 <($hic_read2) \\
            $reads \\
            2> >( tee ${prefix}.stderr.log >&2 )


        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            hifiasm: \$(hifiasm --version 2>&1)
        END_VERSIONS
        """
    } else { // Phasing with Hi-C data is not supported yet
        """
        hifiasm \\
            $args \\
            -o ${prefix}.asm \\
            -t $task.cpus \\
            $reads \\
            2> >( tee ${prefix}.stderr.log >&2 )

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            hifiasm: \$(hifiasm --version 2>&1)
        END_VERSIONS
        """
    }
        stub:
        def args = task.ext.args ?: ''
        def prefix = task.ext.prefix ?: "${meta.id}"
        """
        touch ${prefix}.asm.r_utg.gfa
        touch ${prefix}.asm.ec.bin
        touch ${prefix}.asm.ovlp.source.bin
        touch ${prefix}.asm.ovlp.reverse.bin
        touch ${prefix}.asm.bp.p_ctg.gfa
        touch ${prefix}.asm.p_utg.gfa
        touch ${prefix}.asm.p_ctg.gfa
        touch ${prefix}.asm.a_ctg.gfa
        touch ${prefix}.asm.hap1.p_ctg.gfa
        touch ${prefix}.asm.hap2.p_ctg.gfa
        touch ${prefix}.stderr.log

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            hifiasm: \$(hifiasm --version 2>&1)
        END_VERSIONS
        """

}
