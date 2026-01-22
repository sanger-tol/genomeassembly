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
// scaffolding: scaffolding with YaHS
//
// mito: find the mitogenome using mitohifi
// plastid: find the chloroplast genome using mitohifi
//
// Each stage is defined in the STAGE_SPEC and ORGANELLE_SPEC variables for ease of adding
// new specs later.
//
def stageHifiasmSpec(spec) {
    // Define the datasets types required. We do this additively as once
    // we add a data type, it causes a fork in the specification.
    // bin_data: data required for the production of the overlap graph
    def bin_data = ["long_read_dataset", "long_read_platform"]
    // asm_data: data required for the production of an optionally phased or trio-binned assembly, as well as purging
    def asm_data = bin_data + ["illumina_hic_dataset", "maternal_illumina_dataset", "paternal_illumina_dataset"]
    // polish_scaff_data: data required from scaffolding and beyond
    def polish_scaff_data = asm_data + ["illumina_10x_dataset"]

    // For each stage, define the required datasets, specification parameters, and determine whether a stage
    // should be enabled.
    def STAGE_CONFIG = [
        bin_assembly: [keys: bin_data + ["hifiasm_arguments", "coverage", "busco_lineage"], enabled: true],
        assembly: [keys: asm_data + ["phased_assembly", "trio_assembly", "hifiasm_arguments", "coverage", "busco_lineage"], enabled: true],
        purging: [keys: asm_data + ["purging_cutoffs", "purge_middle", "coverage", "busco_lineage"], enabled: spec.purge],
        polishing: [keys: polish_scaff_data + ["busco_lineage"], enabled: spec.polish],
        scaffolding: [keys: polish_scaff_data + ["yahs_arguments", "busco_lineage"], enabled: spec.scaffold]
    ]

    // Define organelles stages separately as they always depend on the output of the assembly stage
    def ORGANELLE_SPEC = [
        mito: [keys: bin_data + ["assembler", "mitohifi_reference_species", "mitohifi_mito_genetic_code", "mitohifi_arguments"], enabled: spec.find_mito],
        plastid: [keys: bin_data + ["assembler", "mitohifi_reference_species", "mitohifi_plastid_genetic_code", "mitohifi_arguments"], enabled: spec.find_plastid]
    ]

    // Holder variables
    def out_spec = spec.subMap(["id", "assembler"])
    def stageHashes = [:]
    def lastHash = null

    // Process each assembly stage. If the stage is not enabled, we return nothing
    // Otherwise the unique stage hash is calculated using the data, parameters, and if
    // applicable the previous stage's hash.
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

    // Process each organelle stage. If the stage is not enabled, we return nothing
    // Otherwise the hash is calculated using the parameters as well as the stage name,
    // to differentiate organelles and plastids.
    ORGANELLE_SPEC.each { stageName, config ->
        if(!config.enabled) return
        def mode = "contigs"

        // We add the stageName (mito or plastid) to the hash to ensure uniqueness
        def stageSpec = spec.subMap(config.keys)
        def hashContent = stageSpec.values().join("") + stageName + mode
        def prevHash = stageHashes["assembly"]

        def hash = (prevHash + hashContent).sha256()
        stageSpec = stageSpec + [hash: hash, prevHash: prevHash, organelle: stageName, mode: mode]

        out_spec = out_spec + [(stageName): stageSpec]
    }

    return out_spec
}

//
// Process a flat assembly specification for oatk
//
def stageOatkSpec(spec) {
    // Data required for oatk run
    def oatk_data = [
        "long_read_dataset",
        "long_read_platform"
    ]

    // Params required for oatk run
    def oatk_params = [
        "coverage",
        "oatk_kmer_size",
        "oatk_coverage_cutoff",
        "oatk_arguments",
        "oatk_mito_hmm",
        "oatk_plastid_hmm"
    ]

    def oatkSpec = spec.subMap(oatk_data + oatk_params)

    // Generate hash from configuration to detect changes
    def hash_input = oatkSpec.values().join("")
    oatkSpec.hash = hash_input.sha256()

    return spec.subMap(["id", "assembler"]) + [oatk: oatkSpec]
}

//
// Process a flat assembly specification for mitohifi
//
def stageMitohifiSpec(spec) {
    def organelle = "mito"
    def mode = "reads"

    // Data required for mitohifi run
    def mitohifi_data = [
        "long_read_dataset",
        "long_read_platform"
    ]

    // Params required for mitohifi run
    def mitohifi_params = [
        "mitohifi_reference_species",
        "mitohifi_mito_genetic_code",
        "mitohifi_arguments"
    ]

    def mitohifiSpec = spec.subMap(mitohifi_data + mitohifi_params)

    // Generate hash from configuration to detect changes
    def hash_input = mitohifiSpec.values().join("") + organelle + mode
    mitohifiSpec.hash = hash_input.sha256()
    mitohifiSpec.organelle = organelle
    mitohifiSpec.mode = mode

    return spec.subMap(["id", "assembler"]) + [mitohifi: mitohifiSpec]
}
