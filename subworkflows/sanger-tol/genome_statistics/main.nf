include { ASMSTATS            } from '../../../modules/sanger-tol/asmstats/main'
include { BUSCO_BUSCO         } from '../../../modules/nf-core/busco/busco/main'
include { GFASTATS            } from '../../../modules/nf-core/gfastats/main'
include { MERQURYFK_MERQURYFK } from '../../../modules/nf-core/merquryfk/merquryfk/main'

workflow GENOME_STATISTICS {

    take:
    ch_assemblies               // channel: [ val(meta), asm1, asm2 ] - asm2 can be empty
    ch_reads_fastk              // channel: [ val(meta), fastk_hist, [fastk ktabs] ]
    ch_mat_fastk                // channel: [ val(meta), [fastk ktabs] ] - optional
    ch_pat_fastk                // channel: [ val(meta), [fastk ktabs] ] - optional
    val_busco_lineage           // string: busco lineage name
    val_busco_lineage_directory // path: path to local busco lineages directory - optional

    main:
    ch_versions = channel.empty()

    ch_assemblies_split = ch_assemblies
        .flatMap { meta, asm1, asm2 ->
            def meta_asm1 = meta + [_hap: "hap1"]
            def meta_asm2 = meta + [_hap: "hap2"]
            return [ [meta_asm1, asm1], [meta_asm2, asm2] ]
        }
        .filter { _meta, asm -> asm }

    //
    // Module: Calculate assembly stats with asmstats
    //
    ASMSTATS(ch_assemblies_split)
    ch_versions = ch_versions.mix(ASMSTATS.out.versions)

    //
    // Module: Calculate assembly stats with gfastats
    //
    GFASTATS(
        ch_assemblies_split, // assembly
        "fasta",             // out_fmt
        "",                  // genome size
        "",                  // target
        [[],[]],             // agp file
        [[],[]],             // include bed
        [[],[]],             // exclude bed
        [[],[]]              // instructions
    )
    ch_versions = ch_versions.mix(GFASTATS.out.versions)

    //
    // Module: Assess assembly using BUSCO.
    //
    BUSCO_BUSCO(
        ch_assemblies_split,               // assembly
        "genome",                          // busco mode
        val_busco_lineage,                 // lineage to run BUSCO predictions
        val_busco_lineage_directory ?: [], // busco lineage directory
        [],                                // busco config
        true                               // clean intermediates
    )
    ch_versions = ch_versions.mix(BUSCO_BUSCO.out.versions)

    //
    // Module: assess kmer completeness/QV using MerquryFK
    //
    ch_merquryfk_asm_input = ch_assemblies
        .combine(ch_reads_fastk)
        .map { meta_asm, hap1, hap2, _meta_fk, fk_hist, fk_ktabs ->
            [meta_asm, fk_hist, fk_ktabs, hap1, hap2]
        }

    MERQURYFK_MERQURYFK(
        ch_merquryfk_asm_input,
        ch_mat_fastk.ifEmpty([[],[]]),
        ch_pat_fastk.ifEmpty([[],[]])
    )
    ch_versions = ch_versions.mix(MERQURYFK_MERQURYFK.out.versions)

    ch_merquryfk_images = channel.empty()
        .mix(
            MERQURYFK_MERQURYFK.out.spectra_cn_fl,
            MERQURYFK_MERQURYFK.out.spectra_cn_ln,
            MERQURYFK_MERQURYFK.out.spectra_cn_st,
            MERQURYFK_MERQURYFK.out.spectra_asm_fl,
            MERQURYFK_MERQURYFK.out.spectra_asm_ln,
            MERQURYFK_MERQURYFK.out.spectra_asm_st,
            MERQURYFK_MERQURYFK.out.continuity_N,
            MERQURYFK_MERQURYFK.out.block_N,
            MERQURYFK_MERQURYFK.out.block_blob,
            MERQURYFK_MERQURYFK.out.hapmers_blob,
        )
        .transpose()

    emit:
    asmstats             = ASMSTATS.out.stats
    gfastats             = GFASTATS.out.assembly_summary
    busco_summary_txt    = BUSCO_BUSCO.out.short_summaries_txt
    busco_summary_json   = BUSCO_BUSCO.out.short_summaries_json
    merqury_qv           = MERQURYFK_MERQURYFK.out.qv
    merqury_completeness = MERQURYFK_MERQURYFK.out.stats
    merqury_phased_stats = MERQURYFK_MERQURYFK.out.phased_block_stats
    merqury_images       = ch_merquryfk_images
    versions             = ch_versions
}
