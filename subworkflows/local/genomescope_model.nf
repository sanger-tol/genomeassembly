include { CAT_CAT as CAT_CAT_READS   } from "../../modules/nf-core/cat/cat/main"
include { FASTK_FASTK                } from "../../modules/nf-core/fastk/fastk/main"
include { FASTK_HISTEX               } from '../../modules/nf-core/fastk/histex/main'
include { GENESCOPEFK                } from "../../modules/nf-core/genescopefk/main"
include { TRIO_MODE as TRIO_PROCESS  } from '../../subworkflows/local/trio_mode'

workflow GENOMESCOPE_MODEL {

    take:
    reads // [meta, [reads]] 
    matreads // [meta, [matreads]] 
    patreads // [meta, [patreads]]
    trio_flag

    main: 
    ch_versions = Channel.empty()
    matdb_ch = Channel.empty()
    matktab_ch = Channel.empty()
    patdb_ch = Channel.empty()
    patktab_ch = Channel.empty()
    

    trio_flag.view()
    //
    // MODULE: MERGE ALL READS IN ONE FILE
    //
    CAT_CAT_READS( reads )
    ch_versions = ch_versions.mix(CAT_CAT_READS.out.versions)
    
    //
    // LOGIC: KEEP THE CORRECT EXTENSION
    //
    CAT_CAT_READS.out.file_out.map{ meta, reads_ch -> reads_ch.getName().endsWith('gz')
                                    ? [meta, reads_ch.getParent().toString() + '/' + reads_ch.getBaseName().toString() + '.fa.gz'] : [meta, reads_ch.getParent().toString() + '/' + reads_ch.getBaseName().toString() + '.fa'] }
                            .set{ reads_merged_ch }

    //
    // LOGIC: MAKE SURE MERGED READS HAVE THE PROPER EXTENTION
    //
    CAT_CAT_READS.out.file_out.join(reads_merged_ch)
                            .map{ meta, reads_old, reads_new -> reads_old.renameTo(reads_new); }

    //
    // MODULE: GENERATE KMER DATABASE
    //
    FASTK_FASTK( reads_merged_ch )
    ch_versions = ch_versions.mix(FASTK_FASTK.out.versions)

    //
    // MODULE: KEEP THE KMERS HISTOGRAM
    //
    FASTK_HISTEX( FASTK_FASTK.out.hist )
    ch_versions = ch_versions.mix(FASTK_HISTEX.out.versions)


    //
    // MODULE: YAK TO PRODUCE PAT DATABASE
    //

    trio_flag
        .combine( patreads )
        .combine( matreads )
        .branch {
            trio:     it[0] == "trio"
            non_trio: it[0] == "nontrio"
        }
        .set{ trio_data }

    trio_data
            .trio
            .multiMap { trio_mode, pat_meta, pat_data, mat_meta, mat_data ->
                pat:    tuple( pat_meta, pat_data )
                mat:    tuple( mat_meta, mat_data )
            }
            .set{ ch_trio_data }

    TRIO_PROCESS (
        ch_trio_data.pat,
        ch_trio_data.mat
    )    
    
    //
    // MODULE: GENERATE GENOMESCOPE KMER COVERAGE MODEL
    //
    GENESCOPEFK( FASTK_HISTEX.out.hist )
    ch_versions = ch_versions.mix(GENESCOPEFK.out.versions)

    emit:
    model = GENESCOPEFK.out.model
    hist = FASTK_FASTK.out.hist
    ktab = FASTK_FASTK.out.ktab
    pktab = TRIO_PROCESS.out.pktab
    mktab = TRIO_PROCESS.out.mktab
    matdb = TRIO_PROCESS.out.matdb
    patdb = TRIO_PROCESS.out.patdb
    versions = ch_versions
}

