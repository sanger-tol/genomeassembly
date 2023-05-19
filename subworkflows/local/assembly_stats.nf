//
// Evaluate genome by collectings stats
//

//
// Adapted from Sanger Genomenote pipeline by @priyanka-surana
// https://github.com/sanger-tol/genomenote/blob/383f23e6b7a89f9aad6b85c8f7320b5c5825de73/subworkflows/local/genome_statistics.nf
//

include { GFASTATS as GFASTATS_PRI } from '../../modules/nf-core/gfastats/main'
include { GFASTATS as GFASTATS_HAP } from '../../modules/nf-core/gfastats/main'
include { BUSCO               } from '../../modules/nf-core/busco/main'
include { MERQURYFK_MERQURYFK } from '../../modules/nf-core/merquryfk/merquryfk/main'

workflow GENOME_STATISTICS {
    take:
    assembly               // channel: [ meta, primary, haplotigs ]
    lineage                // channel: [ meta, /path/to/buscoDB, lineage ] 
    hist                   // channel: [meta, fastk_hist files]
    ktab                   // channel: [meta, fastk_ktab files]

    main:
    ch_versions = Channel.empty()

    assembly.map{ meta, primary, haplotigs -> [meta, primary] }
        .set{ primary_ch }

    GFASTATS_PRI( primary_ch, 'fasta', [], [], [], [], [], [] )
    ch_versions = ch_versions.mix(GFASTATS_PRI.out.versions.first())
    
    assembly.map{ meta, primary, haplotigs -> [meta, haplotigs] }
        .set{ haplotigs_ch }
    GFASTATS_HAP( haplotigs_ch, 'fasta', [], [], [], [], [], [] )

    // BUSCO
    BUSCO ( assembly.join(lineage)
                    .map{ meta, primary, haplotigs, lineage_db, lineage_name -> 
                            [[id:meta.id, lineage:lineage_name], primary]}, 
            lineage.map{ meta, lineage_db, lineage_name -> lineage_name } ,
            lineage.map{ meta, lineage_db, ch_lineage -> lineage_db },
            [] )
    ch_versions = ch_versions.mix(BUSCO.out.versions.first())
    
    // MerquryFK
    hist.join(ktab).join(assembly)
                    .map{ meta, hist, ktab, primary, hap -> 
                            hap.size() ? [ meta, hist, ktab, primary, hap ] :
                                [ meta, hist, ktab, primary, [] ] } 
                    .set{ ch_merq }
    MERQURYFK_MERQURYFK ( ch_merq )
    ch_versions = ch_versions.mix(MERQURYFK_MERQURYFK.out.versions.first())

    emit:
    busco                     = BUSCO.out.short_summaries_json // meta, path("short_summary.*.json")
    merquryk_completeness     = MERQURYFK_MERQURYFK.out.stats  // meta, stats
    merquryk_qv               = MERQURYFK_MERQURYFK.out.qv     // meta, qv
    assembly_stats_pri        = GFASTATS_PRI.out.assembly_summary  // path("*.assembly_summary")
    assembly_stats_alt        = GFASTATS_HAP.out.assembly_summary  // path("*.assembly_summary")
    versions                  = ch_versions   
}


process GrabFiles {
    tag "${meta.id}"
    executor 'local'

    input:
    tuple val(meta), path("in")

    output:
    tuple val(meta), path("in/*.hist"), path("in/*.ktab*", hidden:true)

    "true"
}
