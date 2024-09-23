include { FASTK_FASTK as FASTK_PAT   } from "../../modules/nf-core/fastk/fastk/main"
include { FASTK_FASTK as FASTK_MAT   } from "../../modules/nf-core/fastk/fastk/main"
include { YAK_COUNT as YAK_COUNT_MAT } from "../../modules/nf-core/yak/count/main"
include { YAK_COUNT as YAK_COUNT_PAT } from "../../modules/nf-core/yak/count/main"

workflow TRIO_MODE {

    take:
    matreads // [meta, [matreads]] 
    patreads // [meta, [patreads]]

    main: 
    ch_versions = Channel.empty()
    
    YAK_COUNT_PAT(patreads)
    FASTK_PAT(patreads)
    patdb_ch = YAK_COUNT_PAT.out.yak
    patktab_ch = FASTK_PAT.out.ktab           
    ch_versions = ch_versions.mix(YAK_COUNT_PAT.out.versions)
    ch_versions = ch_versions.mix(FASTK_PAT.out.versions)
    YAK_COUNT_MAT(matreads)
    FASTK_MAT(matreads)
    matdb_ch = YAK_COUNT_MAT.out.yak
    matktab_ch = FASTK_MAT.out.ktab           
    ch_versions = ch_versions.mix(YAK_COUNT_MAT.out.versions)
    ch_versions = ch_versions.mix(FASTK_MAT.out.versions)

    emit:
    pktab = patktab_ch
    mktab = matktab_ch
    matdb = matdb_ch
    patdb = patdb_ch
    versions = ch_versions
}

