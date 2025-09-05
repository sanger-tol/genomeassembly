include { GAWK as GAWK_GFA_TO_FASTA } from '../../../modules/nf-core/gawk/main'
include { HIFIASM                   } from '../../../modules/nf-core/hifiasm/main'
include { HIFIASM as HIFIASM_BIN    } from '../../../modules/nf-core/hifiasm/main'

workflow RAW_ASSEMBLY {
    take:
    long_reads // channel: [ meta, [reads] ]
    hic_reads  // channel: [ meta, [reads] ]
    trio_dbs   // channel: [ meta, pat_yak, mat_yak ]

    main:
    ch_versions   = Channel.empty()

    //
    // Logic: Hifiasm input channel expects [meta, reads, ul_reads]
    //
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

    //
    // Logic: Get all combinations of input channels:
    //        [long reads, hic/no hic, trio/no trio, bin files]
    //
    //        Then filter out disallowed combinations
    //
    //        Then multimap so that we have the files in the correct channels
    //
    ch_hifiasm_input = long_reads
        | combine(ch_hic_in)
        | combine(ch_trio_in)
        | combine(HIFIASM_BIN.out.bin_files)
        | filter { _lr_meta, _lr, _hic_meta, hic, _trio_meta, pat, mat, _bin_meta, _bin ->
            def is_trio = !(mat.isEmpty() || pat.isEmpty())
            def is_hic  = !hic.isEmpty()
            // Filter disallowed assemblies
            if(is_hic       && is_trio)                     { return false }
            else if(is_hic  && !params.enable_hic_phasing)  { return false }
            else if(is_trio && !params.enable_trio_binning) { return false }
            else                                            { return true  }
        }
        | multiMap { lr_meta, lr, hic_meta, hic, trio_meta, mat, pat, bin_meta, bin ->
            def is_trio = !(mat.isEmpty() || pat.isEmpty())
            def is_hic  = !hic.isEmpty()
            // Add assembly type into the long read meta object
            def assembly_type = "primary"
            if(is_hic)       { assembly_type = "hic_phased"  }
            else if(is_trio) { assembly_type = "trio_binned" }

            def lr_meta_new = lr_meta + [assembly_type: assembly_type, assembly_stage: "raw"]
            long_reads: [lr_meta_new, lr, []]
            hic: [hic_meta, hic]
            trio: [trio_meta, pat, mat]
            bin: [bin_meta, bin]
        }

    //
    // Module: run Hifiasm + resume from bin files
    //
    HIFIASM(
        ch_hifiasm_input.long_reads,
        ch_hifiasm_input.trio,
        ch_hifiasm_input.hic,
        ch_hifiasm_input.bin
    )
    ch_versions = ch_versions.mix(HIFIASM.out.versions)

    //
    // Logic: Mix all the possible Hifiasm GFA assembly channels that
    //        we could be interested in together
    //
    //        Then group them together and filter so that we retain
    //        the correct pair of assemblies in the correct order
    //
    ch_assembly_gfa = Channel.empty()
        | mix(HIFIASM.out.primary_contigs)
        | mix(HIFIASM.out.alternate_contigs)
        | mix(HIFIASM.out.hap1_contigs)
        | mix(HIFIASM.out.hap2_contigs)

    //
    // Module: Convert GFA to FASTA
    //
    GAWK_GFA_TO_FASTA(
        ch_assembly_gfa,
        file("${projectDir}/bin/gfa_to_fasta.awk"),
        false
    )
    ch_versions = ch_versions.mix(GAWK_GFA_TO_FASTA.out.versions)

    //
    // Logic: Split out the correct pri/alt/hap1/hap2 assembly per assembly
    //
    ch_assemblies_gfa = ch_assembly_gfa
        | groupTuple(by: 0, size: 4, remainder: true)
        | map { meta, asms ->
            if(asms.size() == 1) { return null }
            def pri = /hap1.p_ctg.gfa$/
            def alt = /hap2.p_ctg.gfa$/
            if(meta.assembly_type == "primary") {
                if(asms.findAll { it.name =~ /hap/ }.size() == 0) {
                    pri = /^[^.]+\.p_ctg\.gfa$/
                    alt = /a_ctg.gfa$/
                }
            }

            return [meta, asms.find { it.name =~ pri }, asms.find {it.name =~ alt}]
        }
        | filter { it != null }

    ch_assemblies_fasta = GAWK_GFA_TO_FASTA.out.output
        | groupTuple(by: 0, size: 4, remainder: true)
        | map { meta, asms ->
            if(asms.size() == 1) { return null }
            def pri = /hap1.p_ctg.fa$/
            def alt = /hap2.p_ctg.fa$/
            if(meta.assembly_type == "primary") {
                if(asms.findAll { it.name =~ /hap/ }.size() == 0) {
                    pri = /^[^.]+\.p_ctg\.fa$/
                    alt = /a_ctg.fa$/
                }
            }

            return [meta, asms.find { it.name =~ pri }, asms.find {it.name =~ alt}]
        }
        | filter { it != null }

    emit:
    assembly_gfa   = ch_assemblies_gfa
    assembly_fasta = ch_assemblies_fasta
    versions       = ch_versions
}
