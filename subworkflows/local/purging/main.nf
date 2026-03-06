include { FASTA_PURGE_RETAINED_HAPLOTYPE } from '../../../subworkflows/sanger-tol/fasta_purge_retained_haplotype/main'

include { TABIX_BGZIP as BGZIP_PURGED    } from '../../../modules/nf-core/tabix/bgzip/main'

workflow PURGING {
    take:
    ch_purging_specs // channel: spec
    ch_assemblies    // channel: [spec, hap1, hap2]
    val_fastx_reads_per_chunk // int

    main:
    ch_versions = channel.empty()

    //
    // Logic: join all the assemblies with the purging specifications and
    // data, filter for those assemblies which are to be purged, then
    // map out the data for the purging subworkflow.
    //
    ch_purging_input = ch_assemblies
        .combine(ch_purging_specs)
        .filter { asm_meta, _asm1, _asm2, spec -> asm_meta.id == spec.prevID }
        .multiMap { _asm_meta, asm1, asm2, spec ->
            assemblies: [ spec, asm1, asm2 ]
            long_reads: [ spec, spec.data.long_read.reads ]
        }

    //
    // Subworkflow: run purging pipeline
    //
    FASTA_PURGE_RETAINED_HAPLOTYPE(
        ch_purging_input.assemblies,
        ch_purging_input.long_reads,
        val_fastx_reads_per_chunk
    )

    //
    // Module: bgzip all purged assembly fasta
    //
    ch_bgzip_input = FASTA_PURGE_RETAINED_HAPLOTYPE.out.purged_assemblies
        .join(FASTA_PURGE_RETAINED_HAPLOTYPE.out.purged_haplotigs)
        .map { meta, asm1, asm2, haps -> [meta, [asm1, asm2, haps]] }
        .map { meta, fasta ->
            def key = groupKey(meta, fasta.size())
            [ key, fasta ]
        }
        .transpose()

    BGZIP_PURGED(ch_bgzip_input)

    //
    // Logic: combine all purging outputs into a single map for ease of publishing
    //
    ch_purging_output = BGZIP_PURGED.out.output
        .groupTuple(by: 0)
        .map { key, fasta -> [key.target, fasta] }
        .join(FASTA_PURGE_RETAINED_HAPLOTYPE.out.purgedups_splitfa, by: 0)
        .join(FASTA_PURGE_RETAINED_HAPLOTYPE.out.purgedups_splitfa_self_paf, by: 0)
        .join(FASTA_PURGE_RETAINED_HAPLOTYPE.out.purgedups_pbcstat_hist, by: 0)
        .join(FASTA_PURGE_RETAINED_HAPLOTYPE.out.purgedups_pbcstat_basecov, by: 0)
        .join(FASTA_PURGE_RETAINED_HAPLOTYPE.out.purgedups_calcuts_cutoffs, by: 0)
        .join(FASTA_PURGE_RETAINED_HAPLOTYPE.out.purgedups_calcuts_log, by: 0)
        .join(FASTA_PURGE_RETAINED_HAPLOTYPE.out.purgedups_histplot, by: 0)
        .join(FASTA_PURGE_RETAINED_HAPLOTYPE.out.purgedups_bed, by: 0)
        .join(FASTA_PURGE_RETAINED_HAPLOTYPE.out.purgedups_log, by: 0)
        .join(FASTA_PURGE_RETAINED_HAPLOTYPE.out.primary_reads_paf, by: 0)
        .map { spec, fasta, split_fa, splitfa_paf, pbcstat_hist, pbcstat_basecov, calcuts_cutoffs, calcuts_log, histplot, bed, log, reads_paf ->
            return spec.subMap(["id", "stage", "data", "params", "tools"]) + [
                output: [
                    purging: [
                        fasta: fasta,
                        split_fa: split_fa,
                        split_paf: splitfa_paf,
                        pbcstat_hist: pbcstat_hist,
                        pbcstat_basecov: pbcstat_basecov,
                        calcuts_cutoffs: calcuts_cutoffs,
                        calcuts_log: calcuts_log,
                        histplot: histplot,
                        bed: bed,
                        log: log,
                        reads_paf: reads_paf
                    ]
                ]
            ]
        }

    emit:
    purged_assemblies = FASTA_PURGE_RETAINED_HAPLOTYPE.out.purged_assemblies
    purging_output = ch_purging_output
}
