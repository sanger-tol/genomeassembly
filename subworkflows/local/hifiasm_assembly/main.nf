include { HIFIASM                      } from '../../../modules/sanger-tol/hifiasm'
include { HIFIASM as HIFIASM_BIN       } from '../../../modules/sanger-tol/hifiasm'
include { TABIX_BGZIP as BGZIP_HIFIASM } from '../../../modules/nf-core/tabix/bgzip'

workflow HIFIASM_ASSEMBLY {
    take:
    ch_bin_assembly_specs // channel: spec
    ch_assembly_specs     // channel: spec

    main:
    //
    // Module: Run Hifiasm but only generate the .bin files
    //
    ch_hifiasm_bin_input = ch_bin_assembly_specs
        .map { spec ->
            [ spec, spec.data.long_read.reads, [] ]
        }

    HIFIASM_BIN(
        ch_hifiasm_bin_input,
        [[:], []],
        [[:], [], []],
        [[:], []],
        [[:], []]
    )

    //
    // Logic: get the relevant outputs from the Hifiasm bin file generation,
    // then map on the assembly specifications. Filter to match each
    // specification to the required hifiasm outputs, then map out the data.
    //
    ch_hifiasm_input = HIFIASM_BIN.out.bin_files
        .combine(HIFIASM_BIN.out.log, by: 0)
        .combine(ch_assembly_specs)
        .filter { asm_meta, _bin, _log, spec -> asm_meta.id == spec.prevID }
        .multiMap { _asm_meta, bin, log, spec ->
            long_reads: [spec, spec.data.long_read.reads, spec.data.ultralong.reads]
            hic: [spec, spec.params.phased_assembly ? spec.data.hic.reads : []]
            trio: [spec, spec.params.trio_assembly ? spec.data.maternal.yak : [], spec.params.trio_assembly ? spec.data.paternal.yak : []]
            bin: [spec, bin]
            log: [spec, log]
        }

    //
    // Module: run Hifiasm + resume from bin files
    //
    HIFIASM(
        ch_hifiasm_input.long_reads,
        ch_hifiasm_input.hic,
        ch_hifiasm_input.trio,
        ch_hifiasm_input.bin,
        ch_hifiasm_input.log
    )

    //
    // Module: bgzip all hifiasm assembly fasta
    //
    ch_bgzip_input = HIFIASM.out.assembly_fasta
        .map { meta, fasta ->
            def key = groupKey(meta, fasta.size())
            [key, fasta]
        }
        .transpose()

    BGZIP_HIFIASM(ch_bgzip_input)

    //
    // Logic: combine all hifiasm outputs into a single map for ease of publishing
    //
    ch_hifiasm_output = BGZIP_HIFIASM.out.output
        .groupTuple(by: 0)
        .map { key, fasta -> [key.target, fasta] }
        .join(HIFIASM.out.assembly_graphs, by: 0)
        .join(HIFIASM.out.bed, by: 0)
        .join(HIFIASM.out.log, by: 0)
        .map { spec, fasta, graphs, bed, log ->
            return spec.subMap(["id", "stage", "data", "params", "tools"]) + [
                output: [
                    hifiasm: [
                        fasta: fasta,
                        graphs: graphs,
                        bed: bed,
                        log: log
                    ]
                ]
            ]
        }

    //
    // Logic: Split out the correct pri/alt/hap1/hap2 assembly per assembly
    //
    ch_assemblies_fasta = HIFIASM.out.assembly_fasta
        .flatMap { spec, asms ->
            def pri = /hap1.p_ctg.fa$/
            def alt = /hap2.p_ctg.fa$/

            // If hifiasm is run with `--primary` or `-l0`, we do not get
            // hap1/hap2 files - instead we get only p_ctg and a_ctg files, if
            // we don't run in phased or trio mode.
            if(!(asms.find { asm -> asm.name =~ pri } && asms.find { asm -> asm.name =~ alt })) {
                pri = /^[^.]+\.p_ctg\.fa$/
                alt = /a_ctg.fa$/
            }

            return [
                [spec, asms.find { asm -> asm.name =~ pri }, asms.find { asm -> asm.name =~ alt}]
            ]
        }

    emit:
    hifiasm_assemblies = ch_assemblies_fasta
    hifiasm_output = ch_hifiasm_output
}
