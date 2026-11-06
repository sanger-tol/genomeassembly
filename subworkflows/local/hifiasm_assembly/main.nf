include { HIFIASM                            } from '../../../modules/sanger-tol/hifiasm'
include { HIFIASM as HIFIASM_BIN             } from '../../../modules/sanger-tol/hifiasm'
include { HTSLIB_BGZIPTABIX as BGZIP_HIFIASM } from '../../../modules/nf-core/htslib/bgziptabix/main'

workflow HIFIASM_ASSEMBLY {
    take:
    ch_bin_assembly_specs // channel: spec
    ch_assembly_specs // channel: spec

    main:
    //
    // Module: Run Hifiasm but only generate the .bin files
    //
    ch_hifiasm_bin_input = ch_bin_assembly_specs.map { spec ->
        [spec, spec.data.long_read.reads, []]
    }

    HIFIASM_BIN(
        ch_hifiasm_bin_input,
        [[:], []],
        [[:], [], []],
        [[:], []],
        [[:], []],
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
        ch_hifiasm_input.log,
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
        .map { meta, fasta -> [meta, fasta, [], []] }

    BGZIP_HIFIASM(
        ch_bgzip_input,
        "compress",
        false,
        "fa",
    )

    //
    // Logic: combine all hifiasm outputs into a single map for ease of publishing
    //
    ch_hifiasm_output = BGZIP_HIFIASM.out.output
        .groupTuple(by: 0)
        .map { key, fasta -> [key.target, fasta] }
        .join(HIFIASM.out.assembly_graphs, by: 0)
        .join(HIFIASM.out.bed, by: 0)
        .join(HIFIASM.out.log, by: 0)
        .map { spec, fasta, graphs, bed, hifiasm_log ->
            return spec.subMap(["id", "stage", "data", "params", "tools"]) + [
                output: [
                    hifiasm: [
                        fasta: fasta,
                        graphs: graphs,
                        bed: bed,
                        hifiasm_log: hifiasm_log,
                    ],
                ],
            ]
        }

    //
    // Logic: Split out the correct assemblies in sorted order for each Hifiasm run.
    //
    // We prefer hap1/hap2/.../hap_n assemblies if they exist, otherwise we take
    // p_ctg as the primary and a_ctg as the alternative.
    //
    ch_assemblies_fasta = HIFIASM.out.assembly_fasta.map { spec, asms ->
        def validHaps = asms.findAll { asm -> asm =~ /(hap\d+\.p_ctg\.fa$|^.+\.p_ctg\.fa$|^.+\.a_ctg\.fa$)/ }

        // Retain only haplotype assemblies if they exist, otherwise keep all assemblies
        if (validHaps.any { asm -> asm =~ /hap\d+\.p_ctg\.fa$/ }) {
            validHaps = validHaps.findAll { asm -> asm =~ /hap\d+\.p_ctg\.fa$/ }
        }

        def sortedHaps = validHaps.sort { a, b ->
            def extractHapNum = { name ->
                def matcher = name =~ /hap(\d+)\.p_ctg\.fa$/
                matcher.find() ? matcher.group(1).toInteger() : 999999
            }
            def getTigPriority = { name -> name.contains("p_ctg") ? 0 : (name.contains("a_ctg") ? 1 : 2) }

            def numA = extractHapNum.call(a.name)
            def numB = extractHapNum.call(b.name)
            def prioA = getTigPriority.call(a.name)
            def prioB = getTigPriority.call(b.name)

            numA <=> numB ?: prioA <=> prioB ?: a.name <=> b.name
        }

        [spec, sortedHaps]
    }

    emit:
    hifiasm_assemblies = ch_assemblies_fasta
    hifiasm_output     = ch_hifiasm_output
}
