include { ASMSTATS            } from '../../../modules/sanger-tol/asmstats/main'
include { BUSCO_BUSCO         } from '../../../modules/nf-core/busco/busco/main'
include { GFASTATS            } from '../../../modules/nf-core/gfastats/main'
include { MERQURYFK_MERQURYFK } from '../../../modules/nf-core/merquryfk/merquryfk/main'

workflow GENOME_STATISTICS {

    take:
    ch_assemblies               // channel: [ val(meta), asm1, asm2 ] - asm2 can be empty
    ch_fastk                    // channel: [ val(meta), fastk_hist, [fastk ktabs], [mat_fastk_ktabs], [pat_fastk_ktabs] ]
    ch_busco_lineage            // channel: [ val(meta), string: busco_lineage ]
    val_busco_lineage_directory // path: path to local busco lineages directory - optional

    main:
    ch_versions = channel.empty()

    //
    // Logic: rolling check of assembly meta objects to detect duplicates
    //
    def val_asm_meta_list = Collections.synchronizedSet(new HashSet())

    ch_assemblies
        .subscribe { meta, _asm1, _asm2 ->
            if (!val_asm_meta_list.add(meta)) {
                error("Error: Duplicate meta object found in `ch_assemblies` in GENOME_STATISTICS: ${meta}")
            }
        }

    //
    // Logic: split hap1/hap2 into independent channels
    //
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

    //
    // Module: Assess assembly using BUSCO.
    //
    ch_assemblies_for_busco = ch_assemblies
        .map { meta, hap1, hap2 -> [ meta, [hap1, hap2].findAll() ] }
        .join(ch_busco_lineage, by: 0)
        .multiMap { meta, asms, lineage ->
            asms: [ meta, asms ]
            lineage: lineage
        }

    BUSCO_BUSCO(
        ch_assemblies_for_busco.asms,      // assembly
        "genome",                          // busco mode
        ch_assemblies_for_busco.lineage,   // lineage to run BUSCO predictions
        val_busco_lineage_directory ?: [], // busco lineage directory
        [],                                // busco config
        true                               // clean intermediates
    )

    //
    // Module: assess kmer completeness/QV using MerquryFK
    //
    ch_merquryfk_asm_input = ch_assemblies
        .combine(ch_fastk, by: 0)
        .multiMap { meta, hap1, hap2, fk_hist, fk_ktabs, mat_ktabs, pat_ktabs ->
            asms: [meta, fk_hist, fk_ktabs, hap1, hap2]
            mat: [meta, mat_ktabs]
            pat: [meta, pat_ktabs]
        }

    MERQURYFK_MERQURYFK(
        ch_merquryfk_asm_input.asms,
        ch_merquryfk_asm_input.mat,
        ch_merquryfk_asm_input.pat
    )

    //
    // Logic: Join all the outputs into a single map for ease of
    // publishing with workflow outputs
    //
    ch_statistics_output = ASMSTATS.out.stats
        .mix(GFASTATS.out.assembly_summary)
        .map { meta, stats -> [meta - meta.subMap("_hap"), stats] }
        .groupTuple(size: 4)
        .map { meta, out -> [meta, out.flatten().sort { f -> f.getName() }] }

    ch_busco_output = BUSCO_BUSCO.out.batch_summary
        .mix(BUSCO_BUSCO.out.short_summaries_txt)
        .mix(BUSCO_BUSCO.out.short_summaries_json)
        .mix(BUSCO_BUSCO.out.log)
        .mix(BUSCO_BUSCO.out.busco_dir)
        .groupTuple(size: 5)
        .map { meta, out -> [meta, out.flatten().sort { f -> f.getName() }] }

    ch_merqury_output = MERQURYFK_MERQURYFK.out.qv
        .mix(MERQURYFK_MERQURYFK.out.stats)
        .mix(MERQURYFK_MERQURYFK.out.phased_block_stats)
        .mix(MERQURYFK_MERQURYFK.out.images)
        .groupTuple()
        .map { meta, out -> [meta, out.flatten().sort { f -> f.getName() }] }

    emit:
    stats                    = ch_statistics_output
    asmstats                 = ASMSTATS.out.stats
    gfastats                 = GFASTATS.out.assembly_summary
    busco                    = ch_busco_output
    busco_batch_summary      = BUSCO_BUSCO.out.batch_summary
    busco_summary_txt        = BUSCO_BUSCO.out.short_summaries_txt
    busco_summary_json       = BUSCO_BUSCO.out.short_summaries_json
    busco_log                = BUSCO_BUSCO.out.log
    busco_directory          = BUSCO_BUSCO.out.busco_dir
    merqury                  = ch_merqury_output
    merqury_qv               = MERQURYFK_MERQURYFK.out.qv
    merqury_completeness     = MERQURYFK_MERQURYFK.out.stats
    merqury_phased_stats     = MERQURYFK_MERQURYFK.out.phased_block_stats
    merqury_images           = MERQURYFK_MERQURYFK.out.images
}
