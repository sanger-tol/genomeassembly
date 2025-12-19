include { FASTK_HISTEX         } from '../../../modules/nf-core/fastk/histex/main'
include { GENOMESCOPE2         } from "../../../modules/nf-core/genomescope2/main"

include { BUILD_KMER_DATABASES } from '../../../subworkflows/local/build_kmer_databases'

workflow PREPARE_DATA {
    take: 
    ch_specs
    ch_long_reads
    ch_trio
    val_kmer_size

    main:
    ch_versions = channel.empty()

    BUILD_KMER_DATABASES(
        ch_long_reads,
        ch_trio,
        val_kmer_size
    )
    ch_versions = ch_versions.mix(BUILD_KMER_DATABASES.out.versions)

    //
    // Module: FastK histogram to ASCII for Genomescope.
    //         Currently runs only for long reads
    //
    ch_fastk_histex_input = BUILD_KMER_DATABASES.out.ch_long_reads
        .map { meta, _reads, hist, _ktab -> [meta, hist] }

    FASTK_HISTEX(ch_fastk_histex_input)
    ch_versions = ch_versions.mix(FASTK_HISTEX.out.versions)

    //
    // Module: Estimate nuclear coverage with Genomescope
    //
    GENOMESCOPE2(FASTK_HISTEX.out.hist)
    ch_versions = ch_versions.mix(GENOMESCOPE2.out.versions)

    ch_coverage = GENOMESCOPE2.out.model
        .map { meta, model ->
            def kcov_line = model.readLines().find { line -> line =~ /^kmercov/ }
            def kcov = kcov_line ? kcov_line.split(/\s+/).getAt(1).toFloat() : -1
            return [meta, kcov]
        }

    ch_spec_out = ch_spec
        .combine(ch_coverage)
        // .filter { spec, cov_meta, _cov -> spec.long_read_dataset == cov_meta.id && spec.long_read_platform == cov_meta. }
        .map { lr_meta, reads, _hist, _ktab, _cov_meta, cov ->
            def outcov = lr_meta.coverage != -1 ? lr_meta.coverage : cov
            if(outcov == -1) {
                log.error("Error: Unable to get the coverage (either it was not provided in the samplesheet or Genomescope2 failed!")
            }
            def meta_new = lr_meta + [coverage: outcov]
            [meta_new, reads]
        }

    emit:
    ch_specs
    ch_long_reads
    ch_trio
}