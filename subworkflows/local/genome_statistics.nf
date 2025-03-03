//
// Evaluate genome by collectings stats
//

//
// Adapted from Sanger Genomenote pipeline by @priyanka-surana
// https://github.com/sanger-tol/genomenote/blob/383f23e6b7a89f9aad6b85c8f7320b5c5825de73/subworkflows/local/genome_statistics.nf
//

include { GFASTATS as GFASTATS_PRI  } from '../../modules/nf-core/gfastats/main'
include { GFASTATS as GFASTATS_HAP  } from '../../modules/nf-core/gfastats/main'
include { BUSCO as BUSCO_PRI        } from '../../modules/nf-core/busco/main'
include { BUSCO as BUSCO_HAP        } from '../../modules/nf-core/busco/main'
include { MERQURYFK_MERQURYFK       } from '../../modules/nf-core/merquryfk/merquryfk/main'

workflow GENOME_STATISTICS {
    take:
    assembly                  // channel: [ meta, primary, haplotigs ]
    lineage                   // channel: [ meta, /path/to/buscoDB, lineage ] 
    hist                      // channel: [meta, fastk_hist files]
    ktab                      // channel: [fastk_ktab files]
    fastk_pktab               // channel: [fk_pat_ktab files] it is one argument for merquryFK trio case
    phapktab                  // channel: [hapmer files]
    fastk_mktab               // channel: [fk_mat_ktab files] it is one argument for merquryFK trio case
    mhapktab                  // channel: [hapmer files]
    busco_alt                 // channel: true/false

    main:
    ch_versions = Channel.empty()

    //
    // LOGIC: SEPARATE PRIMARY INTO A CHANNEL
    //
    assembly.map{ meta, primary, haplotigs -> [meta, primary] }
        .set{ primary_ch }

    //
    // MODULE: RUN GFASTATS ON PRIMARY ASSEMBLY
    //
    GFASTATS_PRI( primary_ch, 'fasta', [], [], [], [], [], [] )
    ch_versions = ch_versions.mix(GFASTATS_PRI.out.versions.first())
    
    //
    // LOGIC: SEPARATE HAP INTO A CHANNEL
    //
    assembly.map{ meta, primary, haplotigs -> [meta, haplotigs] }
        .set{ haplotigs_ch }

    //
    // MODULE: RUN GFASTATS ON HAPLOTIGS
    //
    GFASTATS_HAP( haplotigs_ch, 'fasta', [], [], [], [], [], [] )

    //
    // MODULE: RUN BUSCO ON PRIMARY ASSEMBLY
    //
    BUSCO_PRI ( primary_ch.join(lineage)
                    .map{ meta, primary, lineage_db, lineage_name -> 
                            [[id:meta.id, lineage:lineage_name], primary]}, 
            lineage.map{ meta, lineage_db, lineage_name -> lineage_name } ,
            lineage.map{ meta, lineage_db, ch_lineage -> lineage_db },
            [] )
    ch_versions = ch_versions.mix(BUSCO_PRI.out.versions.first())
    
    //
    // MODULE: run BUSCO for haplotigs
    // USED FOR HAP1/HAP2 ASSEMBLIES
    //
    if ( busco_alt ) {
        BUSCO_HAP ( haplotigs_ch.join(lineage)
                        .map{ meta, haps, lineage_db, lineage_name -> 
                                [[id:meta.id, lineage:lineage_name], haps]}, 
                lineage.map{ meta, lineage_db, lineage_name -> lineage_name } ,
                lineage.map{ meta, lineage_db, ch_lineage -> lineage_db },
                [] )
    }

    if (params.hifiasm_trio_on) {

        hist.join(ktab).join(assembly)
                       .set{ ch_merq }
        MERQURYFK_MERQURYFK ( ch_merq, fastk_pktab, phapktab, fastk_mktab, mhapktab )

    }
    else{
        //
        // LOGIC: JOIN ASSEMBLY AND KMER DATABASE INPUT
        //
        hist.join(ktab).join(assembly)
                        .map{ meta, hist, ktab, primary, hap -> 
                                hap.size() ? [ meta, hist, ktab, primary, hap ] :
                                    [ meta, hist, ktab, primary, [] ] } 
                        .set{ ch_merq }
        
        //
        // MODULE: RUN KMER ANALYSIS WITH MERQURYFK
        //
        MERQURYFK_MERQURYFK ( ch_merq, [], [], [], [])
    }
    ch_versions = ch_versions.mix(MERQURYFK_MERQURYFK.out.versions.first())

    emit:
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
