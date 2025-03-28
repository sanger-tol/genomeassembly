include { FASTK_FASTK                } from "../../modules/nf-core/fastk/fastk/main"
include { FASTK_HISTEX               } from '../../modules/nf-core/fastk/histex/main'
include { GENESCOPEFK                } from "../../modules/nf-core/genescopefk/main"
include { FASTK_FASTK as FASTK_PAT   } from "../../modules/nf-core/fastk/fastk/main"
include { FASTK_FASTK as FASTK_MAT   } from "../../modules/nf-core/fastk/fastk/main"
include { YAK_COUNT as YAK_COUNT_MAT } from "../../modules/nf-core/yak/count/main"
include { YAK_COUNT as YAK_COUNT_PAT } from "../../modules/nf-core/yak/count/main"
include { MERQURYFK_HAPMAKER         } from "../../modules/local/merquryfk_hapmaker"

workflow KMERS {
    take:
    reads          // [meta, [reads]   ]
    maternal_reads // [meta, [matreads]]
    paternal_reads // [meta, [patreads]]

    main:
    ch_versions = Channel.empty()

    //
    // MODULE: GENERATE KMER DATABASE
    //
    FASTK_FASTK(reads)
    ch_versions = ch_versions.mix(FASTK_FASTK.out.versions)

    //
    // MODULE: KEEP THE KMERS HISTOGRAM
    //
    FASTK_HISTEX( FASTK_FASTK.out.hist )
    ch_versions = ch_versions.mix(FASTK_HISTEX.out.versions)

    //
    // MODULE: GENERATE GENOMESCOPE KMER COVERAGE MODEL
    //
    GENESCOPEFK( FASTK_HISTEX.out.hist )
    ch_versions = ch_versions.mix(GENESCOPEFK.out.versions)

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
    model = GENESCOPEFK.out.model
    hist = FASTK_FASTK.out.hist
    ktab = FASTK_FASTK.out.ktab
    phapktab = TRIO_PROCESS.out.phapktab.ifEmpty( [] )
    mhapktab = TRIO_PROCESS.out.mhapktab.ifEmpty( [] )
    matdb = TRIO_PROCESS.out.matdb.ifEmpty( [] )
    patdb = TRIO_PROCESS.out.patdb.ifEmpty( [] )
    versions = ch_versions
}
