include { FASTK_HISTEX         } from '../../../modules/nf-core/fastk/histex/main'
include { GENOMESCOPE2         } from "../../../modules/nf-core/genomescope2/main"

include { BUILD_KMER_DATABASES } from '../../../subworkflows/local/build_kmer_databases'

workflow PREPARE_DATA {
    take:
    ch_specs
    ch_long_reads
    ch_illumina
    val_kmer_size

    main:
    ch_versions = channel.empty()

    BUILD_KMER_DATABASES(
        ch_specs,
        ch_long_reads,
        ch_illumina,
        val_kmer_size
    )
    ch_versions = ch_versions.mix(BUILD_KMER_DATABASES.out.versions)

    //
    // Module: FastK histogram to ASCII for Genomescope.
    //
    ch_fastk_histex_input = BUILD_KMER_DATABASES.out.long_reads
        .map { meta, _reads, fastk -> [meta, fastk[0]] }

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
            def kcov = kcov_line ? kcov_line.split(/\s+/).getAt(1).toFloat() : null
            return [meta, kcov]
        }

    ch_spec_out = ch_specs
        .combine(ch_coverage)
        .filter { spec, cov_meta, _cov ->
            cov_meta.id == spec.long_read_dataset
        }
        .map { spec, _cov_meta, cov ->
            if(!(spec.long_read_1n_coverage || cov)){
                error(
                    "Error: For assembly ${spec.id}, no long_read_1n_coverage was provided and Genomescope2 failed to " +
                    "estimate the coverage. Please either provide a coverage value, or modify the Genomescope2 parameters " +
                    "and retry."
                )
            }
            return spec + [coverage: spec.long_read_1n_coverage ? spec.long_read_1n_coverage : cov]
        }

    emit:
    specs = ch_spec_out
    long_reads = BUILD_KMER_DATABASES.out.long_reads
    illumina_yakdb = BUILD_KMER_DATABASES.out.illumina_yakdb
    merqury_trio_haptabs = BUILD_KMER_DATABASES.out.merqury_trio_haptabs
    versions = ch_versions
}
