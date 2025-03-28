include { FASTK_FASTK                } from "../../modules/nf-core/fastk/fastk/main"
include { FASTK_HISTEX               } from '../../modules/nf-core/fastk/histex/main'
include { GENESCOPEFK                } from "../../modules/nf-core/genescopefk/main"

workflow KMERS {
    take:
    reads          // [meta, [reads]]
    maternal_reads // [meta, [matreads]]
    paternal_reads // [meta, [patreads]]

    main:
    ch_versions = Channel.empty()

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
    // LOGIC: RUN TRIO WHEN TRIO DATA IS AVAILABLE
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
    //
    // SUBWORKFLOW: RUN TRIO PROCESS WITH TRIO DATA
    //
    TRIO_PROCESS (
        FASTK_FASTK.out.ktab,
        ch_trio_data.mat,
        ch_trio_data.pat
    )
    ch_versions = ch_versions.mix(TRIO_PROCESS.out.versions)

    //
    // MODULE: GENERATE GENOMESCOPE KMER COVERAGE MODEL
    //
    GENESCOPEFK( FASTK_HISTEX.out.hist )
    ch_versions = ch_versions.mix(GENESCOPEFK.out.versions)

    emit:
    model = GENESCOPEFK.out.model
    hist = FASTK_FASTK.out.hist
    ktab = FASTK_FASTK.out.ktab
    phapktab = TRIO_PROCESS.out.phapktab.ifEmpty( [] )
    mhapktab = TRIO_PROCESS.out.mhapktab.ifEmpty( [] )
    matdb = TRIO_PROCESS.out.matdb.ifEmpty( [] )
    patdb = TRIO_PROCESS.out.patdb.ifEmpty( [] )
    versions = ch_versions
}
