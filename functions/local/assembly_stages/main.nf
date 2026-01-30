
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
     data: <map of datatypes to data locations>
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
def stageSpec(spec, paramsConfig, dataList, haptabList) {
    def stageHashes = [:]
    def prevHash = null
    def stages = [:]

    // Create an overall data map
    def dataMap = generateDataMap(spec, dataList, haptabList)

    // Process each assembly stage. If the stage is not enabled, we return nothing
    // Otherwise the unique stage hash is calculated using the data, parameters, and if
    // applicable the previous stage's hash.
    paramsConfig.each { stageName, config ->
        if(!config.enabled) return

        // Extract data and params into submaps
        def datasetsRequired = config.data
        def paramsSpec = spec.subMap(config.params)

        if(config.extraParams) {
            paramsSpec = paramsSpec + config.extraParams
        }

        // Concatenate the data and params values to hash and generate it
        def hashContent = spec.assembler + datasetsRequired.join("") + paramsSpec.values().join("")

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
            dataList: datasetsRequired,
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

    return spec.subMap(["id", "assembler"]) + [hashes: stageHashes, data: dataMap, stages: stages]
}

//
// Add read file and FASTK dataset data to a spec .data list
//
def generateDataMap(spec, dataList, merquryHaptabs) {

    // Define an empty dataset specification
    def emptyDataset = [id: null, platform: null, reads: [], fk_hist: [], fk_ktab: [], yak: [], haptab: []]
    def allDataTypes = ["long_read", "ultralong", "hic", "polishing", "maternal", "paternal"]
    def outputDataMap = allDataTypes.collectEntries { type -> [(type): emptyDataset.clone()] }

    // Get the data types used
    def usedDataTypes = spec.keySet()
        .findAll { k -> k.endsWith("dataset") && spec[k] }
        .collect { type -> type - ~/_dataset/ }

    // For each used data type, find the dataset in the dataList and overwrite
    usedDataTypes.each { dataType ->
        // Get dataset ID and platform
        def dataID = spec[(dataType + "_dataset")]
        def dataPlatform = spec[(dataType + "_platform")]

        // Find dataset
        def dataSet = dataList.find { data ->
            data.id == dataID && data.platform == dataPlatform
        }

        outputDataMap[dataType] = dataSet.clone()
    }

    // Add maternal and paternal haptabs to the dataset
    if(["maternal", "paternal"] in usedDataTypes && merquryHaptabs) {
        def haptabs = merquryHaptabs.find { data ->
            data.long_read_dataset == stageDataMap.long_read_dataset &&
            data.long_read_platform == stageDataMap.long_read_platform &&
            data.maternal_dataset == stageDataMap.maternal_dataset &&
            data.maternal_platform == stageDataMap.maternal_platform &&
            data.paternal_dataset == stageDataMap.paternal_dataset &&
            data.paternal_platform == stageDataMap.paternal_platform
        }

        outputDataMap["maternal"] = outputDataMap["maternal"] + [haptab: haptabs.mat_haptab]
        outputDataMap["paternal"] = outputDataMap["paternal"] + [haptab: haptabs.pat_haptab]
    }

    return outputDataMap
}

//
// Set up a stage for use in a subworkflow. Takes a whole spec and a
// stage name, and returns the stage map with the required dataset attached.
//
def setupStage(spec, stage) {
    if(!spec.stages?[stage]) {
        return null
    }

    def stageData = spec.data.subMap(spec.stages[stage].dataList)

    return spec.stages[stage] + [data: stageData]
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
def stageHifiasmSpec(spec, dataList, haptabList) {
    def DATA_CONFIG = [
        bin_data: ["long_read"],
    ]

    DATA_CONFIG.asm_data = DATA_CONFIG.bin_data + [
        "ultralong",
        "hic",
        "maternal",
        "paternal"
    ]

    DATA_CONFIG.polish_scaff_data = DATA_CONFIG.asm_data + [
        "polishing"
    ]

    def STAGE_CONFIG = [
        bin_assembly: [
            data: DATA_CONFIG.bin_data,
            params: ["hifiasm_bin_arguments", "coverage", "busco_lineage"],
            enabled: true,
            depends: null,
            extraParams: null
        ],
        assembly: [
            data: DATA_CONFIG.asm_data,
            params: ["phased_assembly", "trio_assembly", "hifiasm_bin_arguments", "hifiasm_arguments", "coverage", "busco_lineage"],
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

    return stageSpec(spec, STAGE_CONFIG, dataList, haptabList)
}

/*
    Stage a specification for the Oatk assembly process
*/
def stageOatkSpec(spec, dataList, haptabList) {
    def OATK_CONFIG = [
        oatk: [
            data: ["long_read"],
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

    stageSpec(spec, OATK_CONFIG, dataList, haptabList)
}

/*
    Stage a specification for the Mitohifi assembly process
*/
def stageMitohifiSpec(spec, dataList, haptabList) {
    def MITOHIFI_CONFIG = [
        mitohifi: [
            data: ["long_read"],
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

    return stageSpec(spec, MITOHIFI_CONFIG, dataList, haptabList)
}
