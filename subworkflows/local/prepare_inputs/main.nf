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
        .combine(BUILD_KMER_DATABASES.out.data.collect().map { datasets -> [datasets] })
        .combine(BUILD_KMER_DATABASES.out.merqury_trio_haptabs.collect().map { datasets -> [datasets] }.ifEmpty([[]]) )
        .map { spec, data, haptabs ->
            return addData(spec, data, haptabs)
        }

    emit:
    specs = ch_out_assembly_specs
    versions = ch_versions
}

/* ====================================================
  Assembly specification processing functions

  For simplicity, we pass in the assembly specification in a flat format. We then need to
  post-process each simplification, identifing unique entries at each step of a specification,
  and branch points, to reduce extra computation.

  We identify each independent stage of an assembly using a sha256 hash.
  ===================================================== */

/*
   Process a flat specification for a an assembly, producing a map that nests each
   stage of an assembly specification. Each stage gets a hash that is built from its
   data, params, and the hash of the preceeding stage (if specified). This means that
   the dependency tree of stages can be tracked using md5 sums.

   Takes a paramsConfig - a map of maps, where each key should be the name of a stage
   and each value is a config describing that stage:

   [
     data: <a list of dataset names present in the flat spec>
     params: <a list of param names present in the flat spec>
     enabled: <an expression evaluating to a boolean that enables this stage>
     depends: <the name of a previous stage that this stage depends on, or null>
     extraParams: <a map of extra parameters to add that are not present in the flat spec>
   ]

   Returns a map describing a nested specification:

   [
     id: <assembly id>,
     hashes: <map of key-pair values of stageName -> hash>
     <spec name 1>: <stage map>
     ...
     <spec name N>: <stage map>
   ]

   Where a stage map has the following structure:

   [
     hash: <md5 hash that identifies the stage>
     prevHash: <md5 hash of the stage that this stage depends on>
     data: <map describing the dataset names required by the stage>
     params: <map containing params to parameterise the stage>
   ]

*/
def stageSpec(spec, paramsConfig) {
    def stageHashes = [:]
    def prevHash = null
    def stages = [:]

    // Process each assembly stage. If the stage is not enabled, we return nothing
    // Otherwise the unique stage hash is calculated using the data, parameters, and if
    // applicable the previous stage's hash.
    paramsConfig.each { stageName, config ->
        if(!config.enabled) return

        // Extract data and params into submaps
        def dataSpec = spec.subMap(config.data)
        def paramsSpec = spec.subMap(config.params)

        if(config.extraParams) {
            paramsSpec = paramsSpec + config.extraParams
        }

        // Concatenate the data and params values to hash and generate it
        def hashContent = dataSpec.values().join("") + paramsSpec.values().join("")

        // If a specific dependency is requested, find its hash
        // Otherwise use the last used hash for generation
        def addHash = ""
        if (config.depends) {
            if(!stageHashes[config.depends]) {
                error("Error processing spec [${spec.id}]: No hash has been generated for stage ${config.depends}.")
            }
            addHash = stageHashes[config.depends]
        } else if (prevHash) {
            addHash = prevHash
        }
        hash = (addHash + hashContent).md5()

        // Store the hash in useful places
        stageHashes[stageName] = hash
        stageSpec = [
            stage: stageName,
            id: hash,
            prevID: prevHash,
            data: dataSpec,
            params: paramsSpec
        ]

        // Add the stage to the stages list
        stages = stages + [(stageName): stageSpec]

        // Finally, store the hash for the next iteration but only
        // if we didn't depend on a specific step
        if(!config.depends) {
            prevHash = hash
        }
    }

    return spec.subMap(["id", "assembler"]) + [hashes: stageHashes, stages: stages]
}

/*
    Stage a specification for the Hifiasm assembly process

    Defines three data config levels:
    bin_data: data required for construction of the base overlap graph
    asm_data: data required to produce a genome assembly and purge it
    polish_scaff_data: data requried from the polishing stage onwards

    And 7 assembly stages:
    bin_assembly: construction of the overlap graph
    assembly: production of the (optionally phased, trio) assembly
    purging: purging the assembly
    polishing: polishing the assembly
    scaffolding: scaffolding the assembly
    mito: find mitochondrial genome contigs in the assembly
    plastid: find plastid genome contigs in the assembly
*/
def stageHifiasmSpec(spec) {
    def DATA_CONFIG = [
        bin_data: ["long_read_dataset", "long_read_platform"],
    ]

    DATA_CONFIG.asm_data = DATA_CONFIG.bin_data + [
        "ultralong_dataset",
        "ultralong_platform",
        "hic_dataset",
        "hic_platform",
        "maternal_dataset",
        "maternal_platform",
        "paternal_dataset",
        "paternal_platform"
    ]

    DATA_CONFIG.polish_scaff_data = DATA_CONFIG.asm_data + [
        "polishing_dataset",
        "polishing_platform"
    ]

    def STAGE_CONFIG = [
        bin_assembly: [
            data: DATA_CONFIG.bin_data,
            params: ["hifiasm_arguments", "coverage", "busco_lineage"],
            enabled: true,
            depends: null,
            extraParams: null
        ],
        assembly: [
            data: DATA_CONFIG.asm_data,
            params: ["phased_assembly", "trio_assembly", "hifiasm_arguments", "coverage", "busco_lineage"],
            enabled: true,
            depends: null,
            extraParams: null
        ],
        purging: [
            data: DATA_CONFIG.asm_data,
            params: ["purging_cutoffs", "purge_middle", "coverage", "busco_lineage"],
            enabled: spec.purge,
            depends: null,
            extraParams: null
        ],
        polishing: [
            data: DATA_CONFIG.polish_scaff_data,
            params: ["busco_lineage"],
            enabled: spec.polish,
            depends: null,
            extraParams: null
        ],
        scaffolding: [
            data: DATA_CONFIG.polish_scaff_data,
            params: ["yahs_arguments", "busco_lineage"],
            enabled: spec.scaffold,
            depends: null,
            extraParams: null
        ],
        mito: [
            data: DATA_CONFIG.bin_data,
            params: ["assembler", "mitohifi_reference_species", "mitohifi_mito_genetic_code", "mitohifi_arguments"],
            enabled: spec.find_mito,
            depends: "assembly",
            extraParams: [mode: "contigs", organelle: "mito"]
        ],
        plastid: [
            data: DATA_CONFIG.bin_data,
            params: ["assembler", "mitohifi_reference_species", "mitohifi_plastid_genetic_code", "mitohifi_arguments"],
            enabled: spec.find_plastid,
            depends: "assembly",
            extraParams: [mode: "contigs", organelle: "plastid"]
        ]
    ]

    return stageSpec(spec, STAGE_CONFIG)
}

/*
    Stage a specification for the Oatk assembly process
*/
def stageOatkSpec(spec) {
    def OATK_CONFIG = [
        oatk: [
            data: ["long_read_dataset", "long_read_platform"],
            params: [
                "coverage",
                "oatk_kmer_size",
                "oatk_coverage_cutoff",
                "oatk_arguments",
                "oatk_mito_hmm",
                "oatk_plastid_hmm",
            ],
            enabled: true,
            depends: null,
            extraParams: null
        ]

    ]

    stageSpec(spec, OATK_CONFIG)
}

/*
    Stage a specification for the Mitohifi assembly process
*/
def stageMitohifiSpec(spec) {
    def MITOHIFI_CONFIG = [
        mitohifi: [
            data: ["long_read_dataset", "long_read_platform"],
            params: [
                "mitohifi_reference_species",
                "mitohifi_mito_genetic_code",
                "mitohifi_arguments"
            ],
            enabled: true,
            depends: null,
            extraParams: [mode: "reads", organelle: "mito"]
        ]
    ]

    return stageSpec(spec, MITOHIFI_CONFIG)
}

//
// Add read file and FASTK dataset data to a spec .data list
//
def addData(spec, data_list, merqury_haptabs) {

    // Define an empty dataset specification
    def emptyDataset = [id: null, platform: null, reads: [], fk_hist: [], fk_ktab: [], yak: [], haptab: []]
    def allDataTypes = ["long_read", "ultralong", "hic", "polishing", "maternal", "paternal"]

    // Iterate over the stages in spec and update the data entry
    def specWithData = spec.stages.collectEntries { stageName, stageSpec ->
        // Create a fresh outputDataMap for each stage
        def outputDataMap = allDataTypes.collectEntries { type -> [(type): emptyDataset.clone()] }
        def stageDataMap = stageSpec.data

        // Get the data types used
        def usedDataTypes = stageDataMap.keySet()
            .findAll { k -> k.endsWith("dataset") && stageDataMap[k] }
            .collect { type -> type - ~/_dataset/ }

        usedDataTypes.each { dataType ->
            // Get dataset ID and platform
            def dataID = stageDataMap[(dataType + "_dataset")]
            def dataPlatform = stageDataMap[(dataType + "_platform")]

            // Find dataset
            def dataSet = data_list.find { data ->
                data.id == dataID && data.platform == dataPlatform
            }

            outputDataMap[dataType] = dataSet.clone()
        }

        // Determine if this is a "trio" assembly
        def trioAssembly = ["maternal_dataset", "maternal_platform", "paternal_dataset", "paternal_platform"]
            .every { dataset -> stageDataMap[dataset] }

        // Add maternal and paternal haptabs to the dataset
        if(trioAssembly && merqury_haptabs) {
            def haptabs = merqury_haptabs.find { data ->
                data.long_read_dataset == stageDataMap.long_read_dataset &&
                data.long_read_platform == stageDataMap.long_read_platform &&
                data.maternal_dataset == stageDataMap.maternal_dataset &&
                data.maternal_platform == stageDataMap.maternal_platform &&
                data.paternal_dataset == stageDataMap.paternal_dataset &&
                data.paternal_platform == stageDataMap.paternal_platform
            }

            outputDataMap["maternal"] = outputDataMap["maternal"] + [haptab: haptabs?.mat_haptab]
            outputDataMap["paternal"] = outputDataMap["paternal"] + [haptab: haptabs?.pat_haptab]
        }

        [(stageName): stageSpec + [data: outputDataMap]]
    }

    return spec + [stages: specWithData]
}
