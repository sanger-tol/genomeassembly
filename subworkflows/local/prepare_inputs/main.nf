include { FASTK_HISTEX         } from '../../../modules/nf-core/fastk/histex/main'
include { GENOMESCOPE2         } from "../../../modules/nf-core/genomescope2/main"

include { BUILD_KMER_DATABASES } from '../../../subworkflows/local/build_kmer_databases'

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
        .map { spec ->
            if(spec.assembler == "hifiasm") {
                return stageHifiasmSpec(spec)
            }
            if(spec.assembler == "oatk") {
                return stageOatkSpec(spec)
            }
            if(spec.assembler == "mitohifi") {
                return stageMitohifiSpec(spec)
            }
        }

    emit:
    specs = ch_out_assembly_specs
    data = BUILD_KMER_DATABASES.out.data.collect().map { datasets -> [datasets] }
    merqury_trio_haptabs = BUILD_KMER_DATABASES.out.merqury_trio_haptabs.collect().map { datasets -> [datasets] }
    versions = ch_versions
}

//
// Assembly specification processing functions
//
// For simplicity, we pass in the assembly specification in a flat format. We then need to
// post-process each simplification, identifing unique entries at each step of a specification,
// and branch points, to reduce extra computation.
//
// We identify each independent stage of an assembly using a sha256 hash.
//

//
// Process a flat specification for a Hifiasm assembly, producing the set of specifications
// of each of the following steps:
//
// bin_assembly: production of Hifiasm bin files using long reads only with hifiasm error correction
// assembly: final raw hifiasm assembly, with optional Hi-C phasing or trio binning
// purging: purging pipeline with purge_dups
// polishing: polishing pipeline with longranger + freebayes
// scaffolding: scaffolding with
//
def stageHifiasmSpec(spec) {
    def bin_data = ["long_read_dataset", "long_read_platform"]
    def asm_data = bin_data + ["illumina_hic_dataset", "maternal_illumina_dataset", "paternal_illumina_dataset"]
    def polish_scaff_data = asm_data + ["illumina_10x_dataset"]

    def STAGE_CONFIG = [
        bin_assembly: [keys: bin_data + ["hifiasm_arguments", "coverage"], enabled: true],
        assembly: [keys: asm_data + ["phased_assembly", "trio_assembly", "hifiasm_arguments", "coverage"], enabled: true],
        purging: [keys: asm_data + ["purging_cutoffs", "purge_middle", "coverage"], enabled: spec.purge],
        polishing: [keys: polish_scaff_data + [], enabled: spec.polish],
        scaffolding: [keys: polish_scaff_data + ["yahs_arguments"], enabled: spec.scaffold]
    ]

    def mito_stage = [
        mito: [keys: bin_data + ["find_mito", "find_plastid", "mitohifi_reference_species", "mitohifi_genetic_code", "mitohifi_arguments"], enabled: spec.find_mito || spec.find_plastid]
    ]

    def out_spec = spec.subMap(["id", "assembler"])
    def stageHashes = [:]
    def lastHash = null

    // Process main pipeline
    STAGE_CONFIG.each { stageName, config ->
        if(!config.enabled) return

        def stageSpec = spec.subMap(config.keys)
        def hashContent = stageSpec.values().join("")
        def prevHash = lastHash

        def hash = hashContent.sha256()
        if(prevHash) {
            hash = (prevHash + hashContent).sha256()
        }

        lastHash = hash
        stageHashes[stageName] = hash
        stageSpec = stageSpec + [hash: hash, prevHash: prevHash]

        out_spec = out_spec + [(stageName): stageSpec]
    }

    // Process mito separately (always depends on bin_assembly, independent of pipeline)
    mito_stage.each { stageName, config ->
        if(!config.enabled) return

        def stageSpec = spec.subMap(config.keys)
        def hashContent = stageSpec.values().join("")
        def prevHash = stageHashes["assembly"]

        def hash = (prevHash + hashContent).sha256()
        stageSpec = stageSpec + [hash: hash, prevHash: prevHash]

        out_spec = out_spec + [(stageName): stageSpec]
    }

    return out_spec
}

//
// Process a flat assembly specification for oatk
//
def stageOatkSpec(spec) {
    def oatk_data = ["long_read_dataset", "long_read_platform"]
    def oatk_params = ["coverage", "oatk_kmer_size", "oatk_coverage_cutoff", "oatk_arguments", "oatk_mito_hmm", "oatk_plastid_hmm"]

    def oatkSpec = spec.subMap(oatk_data + oatk_params)
    def oatkHash = oatkSpec.values().join("").sha256()
    oatkSpec = oatkSpec + [hash: oatkHash]

    return spec.subMap(["id", "assembler"]) + [oatk: oatkSpec]
}

//
// Process a flat assembly specification for oatk
//
def stageMitohifiSpec(spec) {
    def mitohifi_data = ["long_read_dataset", "long_read_platform"]
    def mitohifi_params = ["mitohifi_reference_species", "mitohifi_genetic_code", "mitohifi_arguments"]

    def mitohifiSpec = spec.subMap(mitohifi_data + mitohifi_params)
    def mitohifiHash = mitohifiSpec.values().join("").sha256()
    mitohifiSpec = mitohifiSpec + [hash: mitohifiHash]

    return spec.subMap(["id", "assembler"]) + [mitohifi: mitohifiSpec]
}
