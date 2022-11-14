//
// Evaluate genome by collectings stats
//

//
// Adapted from Sanger Genomenote pipeline by @priyanka-surana
// https://github.com/sanger-tol/genomenote/blob/383f23e6b7a89f9aad6b85c8f7320b5c5825de73/subworkflows/local/genome_statistics.nf
//

include { GFASTATS            } from '../../modules/local/gfastats'
include { BUSCO               } from '../../modules/sanger-tol/nf-core-modules/busco/main'
include { MERQURYFK_MERQURYFK } from '../../modules/local/merquryfk'

workflow GENOME_STATISTICS {
    take:
    assembly               // channel: [ meta, primary, haplotigs ]
    lineage             // channel: [ meta, /path/to/buscoDB, lineage ] 
    kmer                   // channel: [ meta, [ /path/to/kmer/kNN ] ]

    main:
    ch_versions = Channel.empty()

     // Get ODB lineage value
    assembly.map{ meta, primary, haplotigs -> [meta, primary] }
        .set{ primary_ch }

    GFASTATS( primary_ch )
    ch_versions = ch_versions.mix(GFASTATS.out.versions.first())

    // BUSCO
    BUSCO ( assembly.map{ meta, primary, haplotigs -> [meta, primary]}, 
            lineage.map{ meta, lineage_db, ch_lineage -> ch_lineage },
            lineage.map{ meta, lineage_db, ch_lineage -> lineage_db },
            [] )
    ch_versions = ch_versions.mix(BUSCO.out.versions.first())
    
    // MerquryFK
    ch_merq = GrabFiles(kmer).combine(assembly).map { meta, hist, ktab, meta2, primary, haplotigs -> [ meta, hist, ktab, primary, haplotigs ] }
    MERQURYFK_MERQURYFK ( ch_merq )
    ch_versions = ch_versions.mix(MERQURYFK_MERQURYFK.out.versions.first())

    emit:
    busco   = BUSCO.out.short_summaries_json // meta, path("short_summary.*.json")
    merquryk_completeness = MERQURYFK_MERQURYFK.out.stats // meta, stats
    merquryk_qv = MERQURYFK_MERQURYFK.out.qv // meta, qv
    assembly_stats = GFASTATS.out.stats // path("*.stats")
    versions = ch_versions   
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
