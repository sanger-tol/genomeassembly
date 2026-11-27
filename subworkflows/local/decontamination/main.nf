include { FASTA_DECONTAMINATION_SYNCON              } from '../../../subworkflows/sanger-tol/fasta_decontamination_syncon/main'

include { HTSLIB_BGZIPTABIX as BGZIP_DECONTAMINATED } from '../../../modules/nf-core/htslib/bgziptabix/main'

workflow DECONTAMINATION {
    take:
    ch_decontamination_specs // channel: spec
    ch_assemblies // channel: [spec, [haps]]

    main:
    //
    // Logic: join all the assemblies with the purging specifications and
    // data, filter for those assemblies which are to be purged, then
    // map out the data for the purging subworkflow.
    //
    ch_decontamination_input = ch_assemblies
        .combine(ch_decontamination_specs)
        .filter { asm_meta, _asms, spec -> asm_meta.id == spec.prevID }
        .multiMap { _asm_meta, asms, spec ->
            assemblies: [spec, asms]
            long_reads: [spec, spec.data.long_read.reads]
        }

    //
    // Subworkflow: run decontamination with Syncon
    //
    FASTA_DECONTAMINATION_SYNCON(
        ch_decontamination_input.assemblies,
        ch_decontamination_input.long_reads,
    )

    //
    // Module: bgzip all purged assembly fasta
    //
    ch_bgzip_input = FASTA_DECONTAMINATION_SYNCON.out.assemblies
        .map { meta, fastas ->
            def key = groupKey(meta, fastas.size())
            [key, fastas]
        }
        .transpose()
        .map { meta, fasta -> [meta, fasta, [], []] }

    BGZIP_DECONTAMINATED(
        ch_bgzip_input,
        "compress",
        false,
        "fa",
    )

    //
    // Logic: combine all decon outputs into a single map for ease of publishing
    //
    ch_decon_output = BGZIP_DECONTAMINATED.out.output
        .groupTuple(by: 0)
        .map { key, fasta -> [key.target, fasta] }
        .join(FASTA_DECONTAMINATION_SYNCON.out.decon_txt)
        .map { spec, fasta, decon_txt ->
            return spec + [output: [decontamination: [
                fasta: fasta,
                decon: decon_txt,
            ]]]
        }

    emit:
    decontaminated_assemblies = FASTA_DECONTAMINATION_SYNCON.out.assemblies
    decontamination_output    = ch_decon_output
}
