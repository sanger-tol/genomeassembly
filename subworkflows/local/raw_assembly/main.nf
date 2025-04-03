include { GFA_TO_FASTA           } from '../../../modules/local/gfa_to_fasta/main'
include { HIFIASM as HIFIASM_BIN } from '../../../modules/nf-core/hifiasm/main'
include { HIFIASM as HIFIASM_HIC } from '../../../modules/nf-core/hifiasm/main'

workflow RAW_ASSEMBLY {
    take:
    long_reads // channel: [ meta, [reads] ]
    hic_reads  // channel: [ meta, [reads] ]
    trio_dbs   // channel: [ meta, pat_yak, mat_yak ]

    main:
    ch_versions   = Channel.empty()
    ch_assemblies = Channel.empty()

    // Hifiasm input channel expects [meta, reads, ul_reads]
    ch_long_reads_input = long_reads
        | map { meta, reads -> [meta, reads, []] }

    //
    // Module: Run Hifiasm but only generate the .bin files
    //
    HIFIASM_BIN(
        ch_long_reads_input,
        [[:], [], []],
        [[:], []],
        [[:], []]
    )
    ch_versions = ch_versions.mix(HIFIASM_BIN.out.versions)

    //
    // Module: Run Hifiasm using previously generated bin files
    //         and generate pri/alt, hic-phased, and trio binned
    //         assemblies IF these are enabled and the data is
    //         available.
    //

    ch_hic_in = hic_reads
        | mix(Channel.of([[:], []]))

    ch_trio_in = trio_dbs
        | mix(Channel.of([[:], [], []]))

    // Get all combinations of input channels:
    //      [long reads, hic/no hic, trio/no trio, bin files]
    // Then filter out disallowed combinations
    // Then multimap so that we have the files in the correct channels
    ch_hifiasm_input = long_reads
        | combine(ch_hic_in)
        | combine(ch_trio_in)
        | combine(HIFIASM_BIN.out.bin_files)
        | filter { lr_meta, lr, hic_meta, hic, trio_meta, pat, mat, bin_meta, bin ->
            // Filter disallowed assemblies
            if(!hic.isEmpty()       && !(mat.isEmpty() && pat.isEmpty()) { return false }
            else if(!hic.isEmpty()  && !params.enable_hic_phasing)       { return false }
            else if(!trio.isEmpty() && !params.enable_trio_binning)      { return false }
            else                                                         { return true  }
        }
        | multimap { lr_meta, lr, hic_meta, hic, trio_meta, mat, pat, bin_meta, bin ->
            // Add assembly type into the long read meta object
            def assembly_type = "primary"
            if(!hic.isEmpty())       { assembly_type = "hic_phased"  }
            else if(!trio.isEmpty()) { assembly_type = "trio_binned" }
            else                     { assembly_type = "unknown"     }
            def lr_meta_new = lr_meta + [assembly_type: assembly_type]

            long_reads: [lr_meta_new, lr, []]
            hic: [hic_meta, hic]
            trio: [trio_meta, pat, mat]
            bin: [bin_meta, bin]
        }

    HIFIASM(
        ch_hifiasm_input.long_reads,
        ch_hifiasm_input.hic,
        ch_hifiasm_input.trio,
        ch_hifiasm_input.bin
    )
    ch_versions = ch_versions.mix(HIFIASM.out.versions)

    // Mix all the possible Hifiasm GFA assembly channels that
    // we could be interested in together
    // Then group them together and filter so that we retain
    // the correct pair of assemblies in the correct order
    ch_assembly_gfa = Channel.empty()
        | mix(HIFIASM.out.primary_contigs)
        | mix(HIFIASM.out.alternate_contigs)
        | mix(HIFIASM.out.hap1_contigs)
        | mix(HIFIASM.out.hap2_contigs)
        | groupTuple(by: 0, size: 4, remainder: true)
        | map { meta, asms ->
            if(asms.size() == 1) { return null }
            def pri = /hap1.p_ctg.gfa$/
            def alt = /hap2.p_ctg.gfa$/
            if(meta.assembly_type == "primary") {
                if(asms.filter { it.name =~ "hap" }.size() == 0) {
                    pri = /^[^.]+\.p_ctg\.gfa$/
                    alt = /a_ctg.gfa$/
                }
            }

            return [meta, [asms.find { it.name =~ pri }, asms.find {it.name =~ alt}]]
        }
        | filter { it != null }
        | transpose
        | map { meta, asm ->
            if(asm.name =~ /hap1.p_ctg.gfa$/ || asm.name =~ /^[^.]+\.p_ctg\.gfa$/) {
                haplotype = "hap1"
            } else {
                haplotype = "hap2"
            }
            def meta_new = meta + [haplotype: haplotype]
            [meta_new, asm]
        }

    //
    // Module: Convert GFA to FASTA and compres with bgzip
    //
    GFA_TO_FASTA(ch_assembly_gfa)
    ch_versions = ch_versions.mix(GFA_TO_FASTA.out.versions)
    ch_assemblies = ch_assemblies
        | mix(GFA_TO_FASTA.out.fasta)

    emit:
    assembly_gfa   = ch_assembly_gfa
    assembly_fasta = ch_assemblies
    versions       = ch_versions
}
