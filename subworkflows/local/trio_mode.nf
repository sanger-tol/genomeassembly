include { FASTK_FASTK as FASTK_PAT   } from "../../modules/nf-core/fastk/fastk/main"
include { FASTK_FASTK as FASTK_MAT   } from "../../modules/nf-core/fastk/fastk/main"
include { YAK_COUNT as YAK_COUNT_MAT } from "../../modules/nf-core/yak/count/main"
include { YAK_COUNT as YAK_COUNT_PAT } from "../../modules/nf-core/yak/count/main"
include { MERQURYFK_HAPMAKER         } from "../../modules/local/hapmaker"

workflow TRIO_MODE {

    take:
    childktab_ch // [meta, [childreads]]
    matreads  // [meta, [matreads]] 
    patreads  // [meta, [patreads]]

    main: 
    ch_versions = Channel.empty()

    //
    // MODULE: GENERATE TRIO DATABASES AND KTABS FOR BOTH PAT AND MAT
    //
    YAK_COUNT_PAT(patreads)
    patdb_ch = YAK_COUNT_PAT.out.yak
    FASTK_PAT(patreads)
    FASTK_PAT.out.ktab.set{patktab_ch}
    ch_versions = ch_versions.mix(YAK_COUNT_PAT.out.versions)
    ch_versions = ch_versions.mix(FASTK_PAT.out.versions)
    
    YAK_COUNT_MAT(matreads)
    matdb_ch = YAK_COUNT_MAT.out.yak 
    FASTK_MAT(matreads)
    FASTK_MAT.out.ktab.set{matktab_ch}
    ch_versions = ch_versions.mix(YAK_COUNT_MAT.out.versions)
    ch_versions = ch_versions.mix(FASTK_MAT.out.versions)

    MERQURYFK_HAPMAKER( matktab_ch, patktab_ch, childktab_ch ) 

    MERQURYFK_HAPMAKER.out.pathap_ktab
    .combine( patktab_ch )
    .map{ hapmeta, pathapktabs, fastkmeta, patktabs -> [hapmeta, pathapktabs, patktabs] }
    .set{pathap_ch}

    MERQURYFK_HAPMAKER.out.mathap_ktab
    .combine( matktab_ch )
    .map{ hapmeta, mathapktabs, fastkmeta, matktabs -> [hapmeta, mathapktabs, matktabs] }
    .set{mathap_ch}

    emit:
    phapktab = pathap_ch
    mhapktab = mathap_ch
    matdb = matdb_ch
    patdb = patdb_ch
    versions = ch_versions
}

