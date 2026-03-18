/**
 * Assembly Specification Processing Functions
 *
 * This module transforms a flat assembly specification into a structured specification
 * for dependency tracking and deduplication. Each assembly stage is uniquely
 * identified via MD5 hash, enabling:
 *
 * - Detection of duplicate work across assembly specifications
 * - Identification of independent stages
 * - Stage-level caching and reuse
 */

/**
 * Transforms a flat assembly specification into a nested structure with stage-specific
 * configurations and dependency tracking via MD5 hashes.
 *
 * This function processes each assembly stage defined in paramsConfig, generating unique
 * MD5 hashes that incorporate the stage's data requirements, parameters, and optionally
 * the hash of a preceding stage. This hash-based approach enables dependency tracking
 * throughout the assembly pipeline.
 *
 * @param spec A flat map containing the assembly specification with id, assembler,
 *             dataset references, and parameter values
 * @param paramsConfig Configuration map defining each assembly stage. Keys are stage names,
 *     values are configuration maps with the following structure:
 *     [
 *       data: <list of dataset names required by this stage>
 *       params: <list of parameter names to extract from spec>
 *       enabled: <boolean expression to enable/disable this stage>
 *       depends: <optional: name of stage this depends on, or null>
 *       extraParams: <optional: additional parameters not in spec>
 *       tools: <optional: list of tools required by this stage>
 *     ]
 * @param dataList List of available dataset objects to resolve data references
 * @param haptabList List of available trio haplotype tables, with information about the
 *     child, mother and father
 *
 * @return A nested map containing the processed assembly specification:
 *     [
 *       id: <original assembly id from spec>
 *       assembler: <assembler name from spec>
 *       hashes: <map of stageName to its MD5 hash string>
 *       data: <consolidated map of all data locations keyed by type>
 *       params: <list of all params used by the spec>
 *       stages: <map of stageName to stage configuration>
 *       tools: <list of all tools used by the spec>
 *     ]
 *
 *     Each stage configuration in the 'stages' map contains:
 *     [
 *       stage: <name of this stage>
 *       id: <unique MD5 hash identifying this stage configuration>
 *       prevID: <MD5 hash of the dependent stage, or empty string if none>
 *       dataList: <list of dataset names this stage requires>
 *       params: <map of parameter names to values for this stage>
 *       tools: <list of all tools used by the stage?
 *     ]
 */
def stageSpec(spec, paramsConfig, dataList, haptabList) {
    // Setup outputs
    def allParams = [:]
    def allTools = [:]
    def allStages = [:]

    // Store hashes - the previous hash, as well as a map of hashes by stages
    def hashesByStage = [:]
    def prevHash = ""

    // Store data - the previous data, as well as a map of data
    def dataByStage = [:]
    def prevData = []

    // Store params - the previous params, as well as a map of params
    def paramsByStage = [:]
    def prevParams = [:]

    // Store tools - the previous tools, as well as a map of tools
    def toolsByStage = [:]
    def prevTools = [:]

    // Create an overall data map
    def dataMap = generateDataMap(spec, dataList, haptabList)

    // Process each assembly stage. If the stage is not enabled, we return nothing
    // Otherwise the unique stage hash is calculated using the data, parameters, and if
    // applicable the previous stage's hash.
    paramsConfig.each { stageName, config ->
        if(!config.enabled) return

        // First, set the stage depenendencies as the previous stage processed.
        //
        // However, if we depend on a named stage, we override this by loading the
        // data for the named stage.
        def dependHash = prevHash
        def dependData = prevData
        def dependParams = prevParams
        def dependTools = prevTools
        if (config.depends) {
            if(!hashesByStage[config.depends] && !paramsByStage[config.depends] && !toolsByStage[config.depends]) {
                error("Error processing spec [${spec.id}]: No stage named ${config.depends} has been generated.")
            }
            dependHash = hashesByStage[config.depends]
            dependData = dataByStage[config.depends]
            dependParams = paramsByStage[config.depends]
            dependTools = toolsByStage.subMap(config.depends)
        }

        // Extract data and params into submaps
        def stageData = dependData + config.data
        def stageParams = dependParams + spec.subMap(config.params)
        def stageTools = dependTools + [(stageName): config.tools]

        // If provided, append the extraParams to the spec
        if(config.extraParams) {
            stageParams = stageParams + config.extraParams
        }

        // Concatenate the data and params values to hash and generate it
        def hashContent = [
            [spec.assembler] +
            stageData +
            stageParams.collect { k, v -> "${k}=${v}" }
        ].join("&")

        def stageHash = (dependHash + hashContent).md5()

        // Save the hash, params and tools in the
        hashesByStage[stageName] = stageHash
        dataByStage[stageName] = stageData
        paramsByStage[stageName] = stageParams
        toolsByStage[stageName] = config.tools

        // Set up the stage specification
        def stageSpec = [
            stage: stageName,
            id: stageHash,
            prevID: dependHash,
            dataList: stageData,
            params: stageParams,
            tools: stageTools
        ]

        // Add the stage to the stages list
        allStages = allStages + [(stageName): stageSpec]

        // Update the list of all params seen
        allParams = allParams + stageParams
        allTools = allTools + stageTools

        // Finally, store the hash for the next iteration but only
        // if we didn't depend on a specific step
        if(!config.depends) {
            prevData = stageData
            prevHash = stageHash
            prevParams = stageParams
            prevTools = stageTools
        }
    }

    return [
        name: spec.id,
        assembler: spec.assembler,
        hashes: hashesByStage,
        data: dataMap,
        params: allParams,
        stages: allStages,
        tools: allTools
    ]
}

/**
 * Generates a consolidated data map by resolving dataset references from the assembly
 * specification and enriching them with haplotype table information when available.
 *
 * This function creates a standardized map containing all data types used in genome
 * assembly (long_read, ultralong, hic, polishing, maternal, paternal). For each data
 * type referenced in the specification, it looks up the corresponding dataset from
 * dataList and populates the output map. Additionally, if maternal and paternal datasets
 * are used, it attempts to locate and attach the appropriate haplotype tables.
 *
 * @param spec The assembly specification containing dataset references in the form
 *             of <datatype>_dataset and <datatype>_platform keys (e.g., long_read_dataset,
 *             long_read_platform)
 * @param dataList List of available dataset objects, each containing:
 *                 [
 *                   id: <dataset identifier>
 *                   platform: <sequencing platform>
 *                   reads: <list of read files>
 *                   fk_hist: <FASTK histogram files>
 *                   fk_ktab: <FASTK k-mer table files>
 *                   yak: <YAK database files>
 *                   haptab: <haplotype table files>
 *                 ]
 * @param merquryHaptabs Optional list of Merqury haplotype table objects that link
 *                       parental datasets to their corresponding haplotype tables.
 *                       Each object should contain dataset/platform identifiers and
 *                       mat_haptab/pat_haptab file references
 *
 * @return A map with keys for all supported data types (long_read, ultralong, hic,
 *         polishing, maternal, paternal). Each key maps to either:
 *         - A populated dataset object if that data type is used in the spec
 *         - An empty dataset object with null/empty values if not used
 *         Maternal and paternal entries will have their haptab field populated if
 *         matching haplotype tables are found in merquryHaptabs
 */
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

        outputDataMap[dataType] = dataSet
    }

    // Add maternal and paternal haptabs to the dataset
    if(merquryHaptabs && ["maternal", "paternal"] in usedDataTypes) {
        def haptabs = merquryHaptabs.find { data ->
            data.long_read_dataset == outputDataMap.long_read_dataset &&
            data.long_read_platform == outputDataMap.long_read_platform &&
            data.maternal_dataset == outputDataMap.maternal_dataset &&
            data.maternal_platform == outputDataMap.maternal_platform &&
            data.paternal_dataset == outputDataMap.paternal_dataset &&
            data.paternal_platform == outputDataMap.paternal_platform
        }

        outputDataMap["maternal"] = outputDataMap["maternal"] + [haptab: haptabs.mat_haptab]
        outputDataMap["paternal"] = outputDataMap["paternal"] + [haptab: haptabs.pat_haptab]
    }

    return outputDataMap
}

/**
 * Configures and stages an assembly specification for the Hifiasm assembler workflow.
 *
 * This function defines the complete pipeline for Hifiasm-based genome assembly, including
 * data requirements and parameters for each stage. It organizes datasets into three tiers
 * based on when they're needed in the pipeline, and configures seven distinct assembly
 * stages with their dependencies and enablement conditions.
 *
 * Assembly stages (in typical execution order):
 * 1. bin_assembly: Constructs the overlap graph from long reads
 * 2. assembly: Produces the genome assembly (optionally phased or trio-based)
 * 3. purging: Removes duplicates and haplotigs (conditional on spec.purge)
 * 4. polishing: Polishes the assembly with additional reads (conditional on spec.polish)
 * 5. scaffolding: Scaffolds contigs using Hi-C data (conditional on spec.scaffold)
 * 6. mito: Identifies mitochondrial contigs (conditional on spec.find_mito, depends on assembly)
 * 7. plastid: Identifies plastid contigs (conditional on spec.find_plastid, depends on assembly)
 *
 * @param spec The flat assembly specification containing assembler parameters, dataset
 *             references, and boolean flags (purge, polish, scaffold, find_mito, find_plastid)
 * @param dataList List of available dataset objects to resolve data references
 * @param haptabList List of available haplotype table objects for trio assembly
 *
 * @return A nested assembly specification with stage configurations and hashes,
 *         as returned by stageSpec()
 */
def stageHifiasmSpec(spec, dataList, haptabList) {
    def STAGE_CONFIG = [
        base: [
            data: ["long_read"],
            params: ["assembler", "hifiasm_bin_arguments", "coverage", "busco_lineage"],
            enabled: true,
            depends: null,
            extraParams: null,
            tools: []
        ],
        hifiasm_assembly: [
            data: [
                "ultralong",
                "hic",
                "maternal",
                "paternal"
            ],
            params: ["phased_assembly", "trio_assembly", "hifiasm_arguments"],
            enabled: true,
            depends: null,
            extraParams: null,
            tools: ["HIFIASM"]
        ],
        purging: [
            data: [],
            params: ["purging_cutoffs", "purge_middle"],
            enabled: spec.purge,
            depends: null,
            extraParams: null,
            tools: ["FASTXALIGN_MINIMAP2ALIGN", "PURGEDUPS_PURGEDUPS"]
        ],
        polishing: [
            data: ["polishing"],
            params: [],
            enabled: spec.polish,
            depends: null,
            extraParams: null,
            tools: ["LONGRANGER_ALIGN", "FREEBAYES"]
        ],
        scaffolding: [
            data: [],
            params: ["yahs_arguments"],
            enabled: spec.scaffold,
            depends: null,
            extraParams: null,
            tools: ["CRAMALIGN_BWAMEM2ALIGNHIC", "CRAMALIGN_MINIMAP2ALIGNHIC", "YAHS"]
        ],
        mitohifi_mito: [
            data: [],
            params: ["mitohifi_reference_species", "mitohifi_mito_genetic_code", "mitohifi_arguments"],
            enabled: spec.find_mito,
            depends: "hifiasm_assembly",
            extraParams: [mode: "contigs", organelle: "mito"],
            tools: ["MITOHIFI_MITOHIFI"]
        ],
        mitohifi_plastid: [
            data: [],
            params: ["mitohifi_reference_species", "mitohifi_plastid_genetic_code", "mitohifi_arguments"],
            enabled: spec.find_plastid,
            depends: "hifiasm_assembly",
            extraParams: [mode: "contigs", organelle: "plastid"],
            tools: ["MITOHIFI_MITOHIFI"]
        ]
    ]

    return stageSpec(spec, STAGE_CONFIG, dataList, haptabList)
}

/**
 * Configures and stages an assembly specification for the Oatk assembler workflow.
 *
 * This function defines a simplified single-stage pipeline for organelle assembly with oatk.
 *
 * Assembly stage:
 * - oatk: Performs organelle assembly from long reads.
 *
 * @param spec The flat assembly specification containing Oatk parameters, dataset
 *             references, and configuration values for k-mer size, coverage cutoffs,
 *             and organelle HMM profiles
 * @param dataList List of available dataset objects to resolve data references
 * @param haptabList List of available haplotype table objects (not used by Oatk but
 *                   required for stageSpec interface compatibility)
 *
 * @return A nested assembly specification with stage configurations and hashes,
 *         as returned by stageSpec()
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
            extraParams: null,
            tools: ["OATK"]
        ]

    ]

    stageSpec(spec, OATK_CONFIG, dataList, haptabList)
}

/**
 * Configures and stages an assembly specification for the MitoHiFi organelle assembler workflow.
 *
 * This function defines a single-stage pipeline for MitoHiFi-based mitochondrial genome assembly.
 *
 * Assembly stage:
 * - mitohifi: Assembles mitochondrial genome directly from long reads
 *
 * @param spec The flat assembly specification containing MitoHiFi parameters, dataset
 *             references, reference species name, and mitochondrial genetic code
 * @param dataList List of available dataset objects to resolve data references
 * @param haptabList List of available haplotype table objects (not used by MitoHiFi but
 *                   required for stageSpec interface compatibility)
 *
 * @return A nested assembly specification with stage configurations and hashes,
 *         as returned by stageSpec()
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
            extraParams: [mode: "reads", organelle: "mito"],
            tools: ["MITOHIFI_MITOHIFI"]
        ]
    ]

    return stageSpec(spec, MITOHIFI_CONFIG, dataList, haptabList)
}

/**
 * Prepares a specific assembly stage for execution by attaching its required datasets.
 *
 * This function extracts a stage configuration from the processed specification and
 * enriches it with the actual dataset objects needed for that stage. It serves as a
 * bridge between the stage specification (which only lists dataset names) and the
 * execution environment (which needs the actual dataset objects).
 *
 * @param spec The processed assembly specification returned by stageSpec(), containing:
 *             - stages: map of stage configurations
 *             - data: consolidated map of all available datasets
 * @param stage The name of the stage to set up (must match a key in spec.stages)
 *
 * @return A stage configuration map ready for execution:
 *         [
 *           stage: <stage name>
 *           id: <stage MD5 hash>
 *           prevID: <previous stage hash or empty string>
 *           dataList: <list of dataset names>
 *           params: <stage parameters>
 *           data: <map of dataset name to dataset object for this stage>
 *         ]
 *         Returns null if the specified stage does not exist in the spec
 */
def setupStage(spec, stage) {
    if (!spec.stages?.get(stage)) {
        return null
    }

    def stageData = spec.data.subMap(spec.stages[stage].dataList)

    return spec.stages[stage] + [data: stageData, assembler: spec.assembler]
}
