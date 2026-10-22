include { FASTK_HISTEX         } from '../../../modules/nf-core/fastk/histex'
include { GENOMESCOPE2         } from "../../../modules/nf-core/genomescope2"
include { GUNZIP               } from '../../../modules/nf-core/gunzip/main'

include { BUILD_KMER_DATABASES } from '../../../subworkflows/local/build_kmer_databases'

include { stageNuclearSpec     } from '../../../functions/local/assembly_stages'
include { stageOrganellarSpec  } from '../../../functions/local/assembly_stages'

workflow PREPARE_INPUTS {
    take:
    ch_specs
    ch_data
    val_kmer_size

    main:
    //
    // Subworkflow: build kmer databases from inputs
    //
    BUILD_KMER_DATABASES(
        ch_specs,
        ch_data,
        val_kmer_size,
    )

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

    ch_coverage = GENOMESCOPE2.out.json_report.map { meta, json ->
        def output = new groovy.json.JsonSlurper().parse(json)
        return [meta, Math.round(output.kcov) ?: null]
    }

    ch_genomescope2_outputs = GENOMESCOPE2.out.linear_plot_png
        .join(GENOMESCOPE2.out.transformed_linear_plot_png, by: 0)
        .join(GENOMESCOPE2.out.log_plot_png, by: 0)
        .join(GENOMESCOPE2.out.transformed_log_plot_png, by: 0)
        .join(GENOMESCOPE2.out.model, by: 0)
        .join(GENOMESCOPE2.out.summary, by: 0)
        .join(GENOMESCOPE2.out.lookup_table, by: 0, remainder: true)
        .join(GENOMESCOPE2.out.fitted_histogram_png, by: 0, remainder: true)
        .join(GENOMESCOPE2.out.json_report, by: 0, remainder: true)
        .map { meta, linear_plot_png, transformed_linear_plot_png, log_plot_png, transformed_log_plot_png, model, summary, lookup_table, fitted_histogram_png, json_report ->
            meta + [
                linear_plot_png: linear_plot_png,
                transformed_linear_plot_png: transformed_linear_plot_png,
                log_plot_png: log_plot_png,
                transformed_log_plot_png: transformed_log_plot_png,
                model: model,
                summary: summary,
                lookup_table: lookup_table,
                fitted_histogram_png: fitted_histogram_png,
                json_report: json_report,
            ]
        }

    //
    // Logic: attach the Genomescope2-calculated coverage to each spec, then either use
    // the provided coverage if provided, otherwise the estimated coverage
    //
    // Also calculate the correct coverage cutoffs for purge_dups and oatk if not provided
    //
    ch_specs_with_coverage = ch_specs
        .combine(ch_coverage)
        .filter { spec, cov_meta, _cov ->
            cov_meta.id == spec.long_read_dataset && cov_meta.platform == spec.long_read_platform
        }
        .map { spec, _cov_meta, cov ->
            if (!(spec.long_read_1n_coverage || cov)) {
                error(
                    """
                    Error: For assembly ${spec.id}, no long_read_1n_coverage was provided and Genomescope2 failed to estimate the coverage.
                    Please either provide a coverage value via long_read_1n_coverage in the assembly specication file, or modify the
                    Genomescope2 parameters with --genomescope2_args and retry.
                    """.stripIndent()
                )
            }
            def coverage = spec.long_read_1n_coverage ? spec.long_read_1n_coverage : cov

            def oatk_coverage_cutoff = spec.oatk_coverage_cutoff ?: coverage * 5
            def purging_cov = coverage + (coverage / 2)
            def purging_cutoffs = spec.purging_cutoffs ?: ",${purging_cov},${purging_cov * 4}"

            return spec + [
                coverage: coverage,
                oatk_coverage_cutoff: oatk_coverage_cutoff,
                purging_cutoffs: purging_cutoffs,
            ]
        }

    //
    // Logic: Pre-computed assemblies can be provided for nuclear assembly if the assembler is set to "input".
    //
    // For convenience we allow these to be zipped, but we need to ensure they are unzipped internally.
    // The following logic unzips any provided zipped assemblies, being very careful to preserve the input file order.
    //
    ch_specs_precomputed_assemblies_to_prepare = ch_specs_with_coverage
        .filter { spec -> spec.assemblies }
        .flatMap { spec ->
            spec.assemblies
                .withIndex()
                .collect { asm, idx ->
                    [spec + [_hap: idx + 1, n_asm: spec.assemblies.size()], asm]
                }
        }
        .branch { _spec, asm ->
            unzip: asm.getExtension() == "gz"
            skip: true
        }

    //
    // Module: Unzip any pre-computed assemblies that are gzipped.
    //
    GUNZIP(ch_specs_precomputed_assemblies_to_prepare.unzip)

    ch_specs_precomputed_assemblies_prepared = GUNZIP.out.gunzip
        .mix(ch_specs_precomputed_assemblies_to_prepare.skip)
        .map { spec, asm ->
            def out_spec = spec - spec.subMap(["_hap", "n_asm"])
            def groupKey = groupKey(out_spec, spec.n_asm)
            [groupKey, [spec, asm]]
        }
        .groupTuple(by: 0)
        .map { key, spec_list ->
            def out_assemblies = spec_list.sort { a, b -> a[0]._hap <=> b[0]._hap }.collect { _spec, asm -> asm }
            return key.getGroupTarget() + [assemblies: out_assemblies]
        }

    ch_specs_to_process = ch_specs_precomputed_assemblies_prepared.mix(ch_specs_with_coverage.filter { spec -> !spec.assemblies })

    //
    // Logic: Prepare assembly specifications
    //
    // For each assembly stage, we create a dedicated meta map that contains the required parameters
    // and data identifiers. We also generate a sha256 hash to uniquely identify each data-parameter
    // combination. Each hashed stage except the first has a dependency on a previous hash, allowing
    // us to fork when required and reduce extraneous computation.
    //
    ch_out_assembly_specs = ch_specs_to_process
        .combine(BUILD_KMER_DATABASES.out.data.collect().map { datasets -> [datasets] })
        .combine(BUILD_KMER_DATABASES.out.merqury_trio_haptabs.collect().map { datasets -> [datasets] }.ifEmpty([[]]))
        .map { spec, data, haptabs ->
            if (spec.assembly_type == "nuclear") {
                return stageNuclearSpec(spec, data, haptabs)
            }
            if (spec.assembly_type == "organellar") {
                return stageOrganellarSpec(spec, data, haptabs)
            }
        }

    emit:
    specs        = ch_out_assembly_specs
    datasets     = BUILD_KMER_DATABASES.out.data
    genomescope2 = ch_genomescope2_outputs
}
