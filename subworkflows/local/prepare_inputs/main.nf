include { FASTK_HISTEX         } from '../../../modules/nf-core/fastk/histex'
include { GENOMESCOPE2         } from "../../../modules/nf-core/genomescope2"

include { BUILD_KMER_DATABASES } from '../../../subworkflows/local/build_kmer_databases'

include { stageHifiasmSpec     } from '../../../functions/local/assembly_stages'
include { stageOatkSpec        } from '../../../functions/local/assembly_stages'
include { stageMitohifiSpec    } from '../../../functions/local/assembly_stages'

workflow PREPARE_INPUTS {
    take:
    ch_specs
    ch_data
    val_kmer_size

    main:
    ch_versions = channel.empty()

    BUILD_KMER_DATABASES(
        ch_specs,
        ch_data,
        val_kmer_size
    )
    ch_versions = ch_versions.mix(BUILD_KMER_DATABASES.out.versions)

    //
    // Module: FastK histogram to ASCII for Genomescope.
    //
    ch_fastk_histex_input = BUILD_KMER_DATABASES.out.data
        .filter { data -> data.platform in ["pacbio_hifi", "oxford_nanopore"] }
        .map { data -> [data.subMap(["id", "platform", "kmer_size"]), data.fk_hist] }

    FASTK_HISTEX(ch_fastk_histex_input)

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

    //
    // Logic: attach the Genomescope2-calculated coverage to each spec, then either use
    // the provided coverage if provided, otherwise the estimated coverage
    //
    ch_specs_with_coverage = ch_specs
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

    //
    // Logic: Prepare assembly specifications
    //
    // For each assembly stage, we create a dedicated meta map that contains the required parameters
    // and data identifiers. We also generate a sha256 hash to uniquely identify each data-parameter
    // combination. Each hashed stage except the first has a dependency on a previous hash, allowing
    // us to fork when required and reduce extraneous computation.
    //
    ch_out_assembly_specs = ch_specs_with_coverage
        .combine(BUILD_KMER_DATABASES.out.data.collect().map { datasets -> [datasets] })
        .combine(BUILD_KMER_DATABASES.out.merqury_trio_haptabs.collect().map { datasets -> [datasets] }.ifEmpty([[]]) )
        .map { spec, data, haptabs ->
            if(spec.assembler == "hifiasm") {
                return stageHifiasmSpec(spec, data, haptabs)
            }
            if(spec.assembler == "oatk") {
                return stageOatkSpec(spec, data, haptabs)
            }
            if(spec.assembler == "mitohifi") {
                return stageMitohifiSpec(spec, data, haptabs)
            }
        }

    emit:
    specs    = ch_out_assembly_specs
    datasets = BUILD_KMER_DATABASES.out.data
    versions = ch_versions
}
