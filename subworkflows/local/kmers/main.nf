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
        | branch { _meta, _reads, hist, _ktab ->
            has_fastk: !hist.isEmpty()
            no_fastk: true
        }

    ch_fastk_input = ch_fastk_status.no_fastk
        | map { meta, reads, _hist, _ktab -> [meta, reads] }

    ch_fastk_skip = ch_fastk_status.has_fastk
        | map { meta, _reads, hist, ktab -> [meta, hist, ktab] }

    FASTK_FASTK(ch_fastk_input)
    ch_versions = ch_versions.mix(FASTK_FASTK.out.versions)

    ch_fastk = FASTK_FASTK.out.hist
        | combine(FASTK_FASTK.out.ktab, by: 0)
        | map { meta, hist, ktab ->
            def meta_new = meta + [kmer_size: params.kmer_size]
            [meta_new, hist, ktab]
        }
        | mix(ch_fastk_skip)

    //
    // Module: FastK histogram to ASCII for Genomescope.
    //         Currently runs only for long reads
    //
    ch_fastk_histex_input = ch_fastk
        | filter { it[0].read_type == "long" && it[0].coverage == -1 }
        | map { meta, hist, _ktab -> [meta, hist] }

    FASTK_HISTEX(ch_fastk_histex_input)
    ch_versions = ch_versions.mix(FASTK_HISTEX.out.versions)

    //
    // Module: Estimate nuclear coverage with Genomescope
    //
    GENOMESCOPE2(FASTK_HISTEX.out.hist)
    ch_versions = ch_versions.mix(GENOMESCOPE2.out.versions)

    ch_coverage = GENOMESCOPE2.out.model
        | map { meta, model ->
            def kcov_line = model.readLines().find { it =~ /^kmercov/ }
            def kcov = kcov_line ? kcov_line.split(/\s+/).getAt(1).toFloat() : -1
            return [meta, kcov]
        }

    ch_long_reads_out = long_reads
        | combine(ch_coverage.ifEmpty([[:], -1]))
        | map { lr_meta, reads, _hist, _ktab, _cov_meta, cov ->
            def outcov = lr_meta.coverage != -1 ? lr_meta.coverage : cov
            if(outcov == -1) {
                log.error("Error: Unable to get the coverage (either it was not provided in the samplesheet or Genomescope2 failed!")
            }
            def meta_new = lr_meta + [coverage: outcov]
            [meta_new, reads]
        }

    //
    // Module: Generate trio databases for maternal and paternal read sets
    //         for trio assembly with hifiasm
    //
    ch_yak_input = paternal_reads.mix(maternal_reads)
        | map { meta, reads, _hist, _ktab -> [meta, reads] }

    YAK_COUNT(ch_yak_input)
    ch_versions = ch_versions.mix(YAK_COUNT.out.versions)

    ch_trio_yak_dbs = YAK_COUNT.out.yak
        | map { meta, yak ->
            def meta_new = meta - meta.subMap("read_type")
            [meta_new, yak]
        }
        | collect
        | map { meta, yaks ->
            def pat = yaks.find { it.name =~ /pat.yak$/ }
            def mat = yaks.find { it.name =~ /mat.yak$/ }
            [meta, pat, mat]
        }

    //
    // Module: Generate trio fastk databases for maternal and paternal read sets
    //         for QC with Merquryfk
    //

    ch_mat_fk   = ch_fastk | filter { it[0].read_type == "mat"  }
    ch_pat_fk   = ch_fastk | filter { it[0].read_type == "pat"  }
    ch_child_fk = ch_fastk | filter { it[0].read_type == "long" }

    MERQURYFK_HAPMAKER(ch_mat_fk, ch_pat_fk, ch_child_fk)
    ch_versions = ch_versions.mix(MERQURYFK_HAPMAKER.out.versions)

    emit:
    long_reads     = ch_long_reads_out
    fastk          = ch_fastk.filter { it[0].read_type == "long" }
    trio_yakdb     = ch_trio_yak_dbs
    maternal_hapdb = MERQURYFK_HAPMAKER.out.mat_hap_ktab
    paternal_hapdb = MERQURYFK_HAPMAKER.out.pat_hap_ktab
    versions       = ch_versions
}
