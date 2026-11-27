include { SYNCON_BUILD as SYNCON_BUILD_ASSEMBLIES } from '../../../modules/nf-core/syncon/build/main'
include { SYNCON_BUILD as SYNCON_BUILD_READS      } from '../../../modules/nf-core/syncon/build/main'
include { SYNCON_DECON                            } from '../../../modules/nf-core/syncon/decon/main'
include { SYNCON_FETCH                            } from '../../../modules/nf-core/syncon/fetch/main'

workflow FASTA_DECONTAMINATION_SYNCON {
    take:
    ch_assemblies // channel: [ val(meta), [fasta1, fasta2, ..., fasta_n] ]
    ch_long_reads // channel: [ val(meta), [ fastx ] ]

    main:
    //
    // Logic: index assemblies by haplotype
    //
    ch_assemblies_indexed = ch_assemblies.flatMap { meta, asms ->
        asms.withIndex().collect { asm, idx -> [meta + [_hap: idx + 1, _nhaps: asms.size()], asm] }
    }

    //
    // Module: Build syncmer database for the reads
    //
    SYNCON_BUILD_READS(ch_long_reads)

    //
    // Module: Build syncmer database for the assemblies
    //
    SYNCON_BUILD_ASSEMBLIES(ch_assemblies_indexed)

    //
    // Module: Identify contaminant contigs in assemblies using the syncmer databases
    //
    ch_syncon_decon_input = SYNCON_BUILD_ASSEMBLIES.out.onekhash
        .join(SYNCON_BUILD_ASSEMBLIES.out.onesyncseq)
        .combine(SYNCON_BUILD_READS.out.onekhash)
        .filter { meta_asm, _asm_1khash, _asm_1syncseq, meta_reads, _reads_1khash ->
            (meta_asm - meta_asm.subMap(["_hap", "_nhaps"])) == meta_reads
        }
        .map { meta_asm, asm_1khash, asm_1syncseq, _meta_reads, reads_1khash ->
            [meta_asm, asm_1khash, asm_1syncseq, reads_1khash]
        }

    SYNCON_DECON(ch_syncon_decon_input, false)

    //
    // Module: Decontaminate the assembly fastas
    //
    SYNCON_FETCH(ch_assemblies_indexed.join(SYNCON_DECON.out.decon))

    //
    // Logic: group the decontaminated assemblies by sample and haplotype.
    // Remove the temporary meta keys _hap and _nhaps.
    //
    ch_output_assemblies = SYNCON_FETCH.out.decontaminated_fasta
        .map { meta, asm ->
            def out_meta = meta - meta.subMap(["_hap", "_nhaps"])
            def groupKey = groupKey(out_meta, meta._nhaps)
            [groupKey, [meta, asm]]
        }
        .groupTuple(by: 0)
        .map { key, meta_asm_list ->
            def out_assemblies = meta_asm_list.sort { a, b -> a[0]._hap <=> b[0]._hap }.collect { _spec, asm -> asm }
            return [key.getGroupTarget(), out_assemblies]
        }

    ch_output_decon_txt = SYNCON_DECON.out.decon
        .map { meta, decon_txt ->
            def out_meta = meta - meta.subMap(["_hap", "_nhaps"])
            def groupKey = groupKey(out_meta, meta._nhaps)
            [groupKey, [meta, decon_txt]]
        }
        .groupTuple(by: 0)
        .map { key, meta_decon_list ->
            def out_decon_txt = meta_decon_list.sort { a, b -> a[0]._hap <=> b[0]._hap }.collect { _spec, decon_txt -> decon_txt }
            return [key.getGroupTarget(), out_decon_txt]
        }

    emit:
    assemblies = ch_output_assemblies // channel: [ val(meta), [fasta1, fasta2, ..., fasta_n] ]
    decon_txt  = ch_output_decon_txt // channel: [ val(meta), [decon_txt1, decon_txt2, ..., decon_txt_n] ]
}
