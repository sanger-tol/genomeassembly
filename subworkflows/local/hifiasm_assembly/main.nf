include { HIFIASM                      } from '../../../modules/sanger-tol/hifiasm/main'
include { HIFIASM as HIFIASM_BIN       } from '../../../modules/sanger-tol/hifiasm/main'
include { TABIX_BGZIP as BGZIP_HIFIASM } from '../../../modules/nf-core/tabix/bgzip/main'

workflow HIFIASM_ASSEMBLY {
    take:
    ch_bin_assembly_specs // channel: spec
    ch_assembly_specs     // channel: spec
    ch_data               // channel: data_map

    main:
    ch_versions = channel.empty()

    //
    // Module: Run Hifiasm but only generate the .bin files
    //
    ch_hifiasm_bin_input = ch_bin_assembly_specs
        .combine(ch_data)
        .map { spec, datasets ->
            def long_reads = datasets.find { data ->
                data.id == spec.long_read_dataset && data.platform == spec.long_read_platform
            }

            [ spec, long_reads.reads, [] ]
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
    // then map on the assembly specifications and data. Filter to match each
    // specification to the required hifiasm outputs, then map out the data.
    //
    ch_hifiasm_input = HIFIASM_BIN.out.bin_files
        .combine(HIFIASM_BIN.out.log, by: 0)
        .combine(ch_assembly_specs)
        .combine(ch_data)
        .filter { meta, _bin, _log, spec, _datasets -> meta.hash == spec.prevHash }
        .multiMap { meta, bin, log, spec, datasets ->
            def long_read_data = datasets.find { data ->
                data.id == spec.long_read_dataset && data.platform == spec.long_read_platform
            }
            def hic_data = datasets.find { data ->
                data.id == spec.illumina_hic_dataset && data.platform == "illumina_hic"
            }
            def maternal_data = datasets.find { data ->
                data.id == spec.illumina_hic_dataset && data.platform == "illumina"
            }
            def paternal_data = datasets.find { data ->
                data.id == spec.illumina_hic_dataset && data.platform == "illumina"
            }

            def hic_reads = hic_data?.reads && spec.phased_assembly ? hic_data.reads : []
            def mat_yak = maternal_data?.yak && spec.trio_assembly? maternal_data.yak : []
            def pat_yak = paternal_data?.yak && spec.trio_assembly? paternal_data.yak : []

            long_reads: [ spec, long_read_data.reads, [] ]
            hic: [ spec, hic_reads ]
            trio: [ spec, mat_yak, pat_yak ]
            bin: [ spec, bin ]
            log: [ spec, log ]
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
            [ key, fasta ]
        }
        .transpose()

    BGZIP_HIFIASM(ch_bgzip_input)
    ch_versions = ch_versions.mix(BGZIP_HIFIASM.out.versions)

    //
    // Logic: combine all hifiasm outputs into a single map for ease of publishing
    //
    ch_hifiasm_output = BGZIP_HIFIASM.out.output
        .groupTuple(by: 0)
        .map { key, fasta -> [key.target, fasta] }
        .join(HIFIASM.out.assembly_graphs, by: 0)
        .join(HIFIASM.out.bed, by: 0)
        .join(HIFIASM.out.log, by: 0)
        .join(HIFIASM.out.corrected_reads, remainder: true, by: 0)
        .join(HIFIASM.out.read_overlaps, remainder: true, by: 0)
        .map { spec, fasta, graphs, bed, log, reads, paf ->
            return [
                id: spec.hash,
                fasta: fasta,
                graphs: graphs,
                bed: bed,
                log: log,
                reads: reads,
                paf: paf
            ]
        }

    //
    // Logic: Split out the correct pri/alt/hap1/hap2 assembly per assembly
    //
    ch_assemblies_fasta = HIFIASM.out.assembly_fasta
        .flatMap { meta, asms ->
            def pri = /hap1.p_ctg.fa$/
            def alt = /hap2.p_ctg.fa$/

            // If hifiasm is run with `--primary` or `-l0`, we do not get
            // hap1/hap2 files - instead we get only p_ctg and a_ctg files, if
            // we don't run in phased or trio mode.
            def primary_flags = meta.hifiasm_arguments =~ /(--primary)|(-l\s*0)/
            if(primary_flags && !(meta.phased_assembly || meta.trio_assembly) ) {
                pri = /^[^.]+\.p_ctg\.fa$/
                alt = /a_ctg.fa$/
            }

            return [
                [meta, asms.find { asm -> asm.name =~ pri }, asms.find { asm -> asm.name =~ alt}]
            ]
        }

    emit:
    hifiasm_assemblies = ch_assemblies_fasta
    hifiasm_all = ch_hifiasm_output
    versions = ch_versions
}
