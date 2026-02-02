include { GAWK as GAWK_GFA_TO_FASTA } from '../../../modules/nf-core/gawk/main'
include { HIFIASM                   } from '../../../modules/sanger-tol/hifiasm/main'
include { HIFIASM as HIFIASM_BIN          } from '../../../modules/sanger-tol/hifiasm/main'
include { CAT_CAT as CONCATENATE_ALTERNATES } from '../../../modules/nf-core/cat/cat'

workflow RAW_ASSEMBLY {
    take:
    long_reads // channel: [ meta, [reads] ]
    hic_reads  // channel: [ meta, [reads] ]
    trio_dbs   // channel: [ meta, pat_yak, mat_yak ]

    main:
    ch_versions   = channel.empty()

    //
    // Logic: Hifiasm input channel expects [meta, reads, ul_reads]
    //
    ch_long_reads_input = long_reads
        .map { meta, reads -> [meta, reads, []] }

    //
    // Module: Run Hifiasm but only generate the .bin files
    //
    HIFIASM_BIN(
        ch_long_reads_input,
        [[:], []],
        [[:], [], []],
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
        .mix(channel.of([[:], []]))

    ch_trio_in = trio_dbs
        .mix(channel.of([[:], [], []]))

    //
    // Logic: Get all combinations of input channels:
    //        [long reads, hic/no hic, trio/no trio, bin files]
    //
    //        Then filter out disallowed combinations
    //
    //        Then multimap so that we have the files in the correct channels
    //
    ch_hifiasm_input = long_reads
        .combine(ch_hic_in)
        .combine(ch_trio_in)
        .combine(HIFIASM_BIN.out.bin_files)
        .filter { _lr_meta, _lr, _hic_meta, hic, _trio_meta, pat, mat, _bin_meta, _bin ->
            def is_trio = !(mat.isEmpty() || pat.isEmpty())
            def is_hic  = !hic.isEmpty()
            // Filter disallowed assemblies
            if(is_hic       && is_trio)                     { return false }
            else if(is_hic  && !params.enable_hic_phasing)  { return false }
            else if(is_trio && !params.enable_trio_binning) { return false }
            else                                            { return true  }
        }
        .multiMap { lr_meta, lr, hic_meta, hic, trio_meta, mat, pat, bin_meta, bin ->
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
        ch_hifiasm_input.hic,
        ch_hifiasm_input.trio,
        ch_hifiasm_input.bin
    )
    ch_versions = ch_versions.mix(HIFIASM.out.versions)

    //
    // Logic: Split out the correct pri/alt/hap1/hap2 assembly per assembly
    //
    ch_assemblies_fasta = HIFIASM.out.assembly_fasta
        .transpose()
        .groupTuple(by: 0)
        .map { meta, asms ->
            // Identify Primary (hap1.p_ctg.fa or p_ctg.fa)
            def pri_pattern_hap = /hap1.p_ctg.fa$/
            def pri_pattern_gen = /^[^.]+\.p_ctg\.fa$/

            // Filter to only keep contig files (*ctg.fa)
            def ctg_asms = asms.findAll { it.name =~ /ctg.fa/ }

            def pri = ctg_asms.find { it.name =~ pri_pattern_hap }
            if (!pri) {
                 // Fallback: check generally for p_ctg.fa if no "hap" naming is used in ANY file
                 def has_hap = ctg_asms.any { it.name =~ /hap/ }
                 if (!has_hap) {
                      pri = ctg_asms.find { it.name =~ pri_pattern_gen }
                 }
            }

            // Alternates are everything else
            def alts = ctg_asms - pri

            [ meta, pri, alts ]
        }
        .branch { meta, pri, alts ->
             concat: alts.size() > 1
             ready: true
        }

    CONCATENATE_ALTERNATES(
        ch_assemblies_fasta.concat.map { meta, pri, alts -> [ meta, alts ] }
    )
    ch_versions = ch_versions.mix(CONCATENATE_ALTERNATES.out.versions)

    ch_assemblies_ready = ch_assemblies_fasta.concat
        .join(CONCATENATE_ALTERNATES.out.file_out)
        .map { meta, pri, alts_orig, merged_alt -> [ meta, pri, [merged_alt] ] }
        .mix(ch_assemblies_fasta.ready)
        .map { meta, pri, alts ->
            def alt = alts ? alts[0] : []
            [ meta, pri, alt ]
        }

    emit:
    hifiasm_fasta = ch_assemblies_ready
    hifiasm_gfa   = HIFIASM.out.assembly_graphs
    hifiasm_bed   = HIFIASM.out.bed
    hifiasm_log   = HIFIASM.out.log
    versions      = ch_versions
}
