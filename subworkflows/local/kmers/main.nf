include { FASTK_FASTK        } from "../../../modules/nf-core/fastk/fastk/main"
include { FASTK_HISTEX       } from '../../../modules/nf-core/fastk/histex/main'
include { GENOMESCOPE2       } from "../../../modules/nf-core/genomescope2/main"
include { YAK_COUNT          } from "../../../modules/nf-core/yak/count/main"
include { MERQURYFK_HAPMAKER } from "../../../modules/nf-core/merquryfk/hapmaker"

workflow KMERS {
    take:
    long_reads     // [meta, [reads], fk_hist, [fk_ktab]]
    maternal_reads // [meta, [reads], fk_hist, [fk_ktab]]
    paternal_reads // [meta, [reads], fk_hist, [fk_ktab]]

    main:
    ch_versions = Channel.empty()

    //
    // Module: Generate FastK databases for all read sets without one
    //
    ch_fastk_status = long_reads.mix(maternal_reads, paternal_reads)
        | branch { meta, reads, hist, ktab ->
            has_fastk: !hist.isEmpty()
            no_fastk: true
        }

    ch_fastk_status.no_fastk
        | map { meta, reads, _hist, _ktab -> [meta, reads] }
        | FASTK_FASTK
    ch_versions = ch_versions.mix(FASTK_FASTK.out.versions)

    ch_fastk = FASTK_FASTK.out.hist
        | combine(FASTK_FASTK.out.ktab, by: 0)
        | mix(ch_fastk_status.has_fastk)

    //
    // Module: FastK histogram to ASCII for Genomescope
    //
    ch_fastk
        | filter { it[0].reads == "long" }
        | map { meta, hist, _ktab -> [meta, hist] }
        | FASTK_HISTEX
    ch_versions = ch_versions.mix(FASTK_HISTEX.out.versions)

    //
    // Module: Estimate nuclear coverage with Genomescope
    //
    GENOMESCOPE2(FASTK_HISTEX.out.hist)
    ch_versions = ch_versions.mix(GENESCOPEFK.out.versions)

    ch_coverage = GENESCOPEFK.out.model
        | map { meta, model ->
            def kcov_line = model.readLines().find { it =~ "kmercov" }
            def kcov = kcov_line ? kcov_line.tokenize(/\s+/).getAt(1).toFloat() : -1
            return [meta, kcov]
        }

    //
    // Module: Generate trio databases for maternal and paternal read sets
    //         for trio assembly with hifiasm
    //
    reads
        | filter { it[0].reads in ["pat", "mat"] }
        | map { meta, reads, _hist, _ktab -> [meta, reads] }
        | YAK_COUNT
    ch_versions = ch_versions.mix(YAK_COUNT.out.versions)

    ch_trio_yak_dbs = YAK_COUNT.out.yak
        | map { meta, yak ->
            def meta_new = meta - meta.subMap("trio")
            [meta_new, yak]
        }
        | collect()
        | map { meta, yaks ->
            def pat = yaks.find { it.name =~ "pat" }
            def mat = yaks.find { it.name =~ "mat" }
            [meta, pat, mat]
        }

    //
    // Module: Generate trio fastk databases for maternal and paternal read sets
    //         for QC with Merquryfk
    //

    ch_mat_fk   = ch_fastk | filter { it[0].reads == "mat"  }
    ch_pat_fk   = ch_fastk | filter { it[0].reads == "pat"  }
    ch_child_fk = ch_fastk | filter { it[0].reads == "long" }

    MERQURYFK_HAPMAKER(ch_mat_fk, ch_pat_fk, ch_child_fk)
    ch_versions = ch_versions.mix(MERQURYFK_HAPMAKER.out.versions)

    ch_trio_hap_dbs = MERQURYFK_HAPMAKER.out.pat_hap_ktab
        | combine(MERQURYFK_HAPMAKER.out.mat_hap_ktab, by: 0)

    emit:
    coverage   = ch_coverage
    fastk      = ch_fastk.filter { it[0].reads == "long" }
    trio_yakdb = ch_trio_yak_dbs
    trio_hapdb = ch_trio_hap_dbs
    versions   = ch_versions
}
