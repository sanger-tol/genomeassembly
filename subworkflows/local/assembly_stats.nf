//
// Evaluate genome by collectings stats
//

//
// Adapted from Sanger Genomenote pipeline by @priyanka-surana
// https://github.com/sanger-tol/genomenote/blob/383f23e6b7a89f9aad6b85c8f7320b5c5825de73/subworkflows/local/genome_statistics.nf
//

include { BUSCO    } from '../../modules/sanger-tol/nf-core-modules/busco/main'
include { MERQURYFK_MERQURYFK } from '../../modules/local/merquryfk'
include { GET_ODB             } from '../../modules/local/get_odb'

workflow GENOME_STATISTICS {
    take:
    assembly               // channel: [ meta, primary, haplotigs ]
    lineage_db             // channel: /path/to/buscoDB
    kmer                   // channel: [ meta, [ /path/to/kmer/kNN ] ]

    main:
    ch_versions = Channel.empty()

     // Get ODB lineage value
    assembly.map{ meta, primary, haplotigs -> [meta, primary] }
        .set{ primary_ch }
    GET_ODB ( primary_ch )
    ch_versions = ch_versions.mix(GET_ODB.out.versions.first())

    // BUSCO
    ch_lineage = GET_ODB.out.csv.splitCsv().map { row -> row[1] }
    BUSCO ( assembly.map{ meta, primary, haplotigs -> [meta, primary]}, ch_lineage, lineage_db, [] )
    ch_versions = ch_versions.mix(BUSCO.out.versions.first())
    
    // MerquryFK
    // TODO make merquryfk accept haplotigs
    ch_merq = GrabFiles(kmer).combine(assembly).map { meta, hist, ktab, meta2, primary, haplotigs -> [ meta, hist, ktab, primary, haplotigs ] }
    MERQURYFK_MERQURYFK ( ch_merq )
    ch_versions = ch_versions.mix(MERQURYFK_MERQURYFK.out.versions.first())

    emit:
    BUSCO   = BUSCO.out.short_summaries_json // meta, path("short_summary.*.json")
    MERQURY_COMPLETENESS = MERQURYFK_MERQURYFK.out.stats // meta, stats
    MERQURY_QV = MERQURYFK_MERQURYFK.out.qv // meta, qv

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
