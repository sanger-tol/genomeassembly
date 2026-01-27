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
   Process a flat specification for a Hifiasm assembly, producing the set of specifications
   of each of the following steps:

   bin_assembly: production of Hifiasm bin files using long reads only with hifiasm error correction
   assembly: final raw hifiasm assembly, with optional Hi-C phasing or trio binning
   purging: purging pipeline with purge_dups
   polishing: polishing pipeline with longranger + freebayes
   scaffolding: scaffolding with YaHS

   mito: find the mitogenome using mitohifi
   plastid: find the chloroplast genome using mitohifi

   Each stage is defined in the STAGE_SPEC and ORGANELLE_SPEC variables for ease of adding
   new specs later.
*/
def stageHifiasmSpec(spec) {
    // Define the datasets types required. We do this additively as once we add a data type,
    // it causes a fork in the specification.
    //
    // bin_data: data required for the production of the overlap graph
    // asm_data: data required for assembly and purging
    // polish_scaff_data: data required for purging and scaffolding
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

    // Define the data and parameters needed at each step of the nuclear assembly process.
    //
    // We define each stage as a map with a list containing the required fields for the data,
    // a list of the required fields for parameterisation, and a flag that determines whether
    // the stage is enabled or not.
    def STAGE_CONFIG = [
        bin_assembly: [
            data: DATA_CONFIG.bin_data,
            params: ["hifiasm_arguments", "coverage", "busco_lineage"],
            enabled: true
        ],
        assembly: [
            data: DATA_CONFIG.asm_data,
            params: ["phased_assembly", "trio_assembly", "hifiasm_arguments", "coverage", "busco_lineage"],
            enabled: true
        ],
        purging: [
            data: DATA_CONFIG.asm_data,
            params: ["purging_cutoffs", "purge_middle", "coverage", "busco_lineage"],
            enabled: spec.purge
        ],
        polishing: [
            data: DATA_CONFIG.polish_scaff_data,
            params: ["busco_lineage"],
            enabled: spec.polish
        ],
        scaffolding: [
            data: DATA_CONFIG.polish_scaff_data,
            params: ["yahs_arguments", "busco_lineage"],
            enabled: spec.scaffold
        ]
    ]

    // Define the data and parameters needed at each step of the organellar  assembly process.
    //
    // We define each stage as a map with a list containing the required fields for the data,
    // a list of the required fields for parameterisation, and a flag that determines whether
    // the stage is enabled or not.
    def ORGANELLE_SPEC = [
        mito: [
            data: DATA_CONFIG.bin_data,
            params: ["assembler", "mitohifi_reference_species", "mitohifi_mito_genetic_code", "mitohifi_arguments"],
            enabled: spec.find_mito
        ],
        plastid: [
            data: DATA_CONFIG.bin_data,
            params: ["assembler", "mitohifi_reference_species", "mitohifi_plastid_genetic_code", "mitohifi_arguments"],
            enabled: spec.find_plastid
        ]
    ]

    // Holder variables
    def output_specification = spec.subMap(["id", "assembler"])
    def stageHashes = [:]
    def prevHash = null

    // Process each assembly stage. If the stage is not enabled, we return nothing
    // Otherwise the unique stage hash is calculated using the data, parameters, and if
    // applicable the previous stage's hash.
    STAGE_CONFIG.each { stageName, config ->
        if(!config.enabled) return

        // Extract data and params into submaps
        def dataSpec = spec.subMap(config.data)
        def paramsSpec = spec.subMap(config.params)

        // Concatenate the data and params values to hash and generate it
        // If we have a previous hash, we add this to the computed hash and recompute
        def hashContent = dataSpec.values().join("") + paramsSpec.values().join("")
        def hash = hashContent.sha256()

        if(prevHash) {
            hash = (prevHash + hashContent).sha256()
        }

        // Store the hash in useful places
        stageHashes[stageName] = hash
        stageSpec = [
            hash: hash,
            prevHash: prevHash,
            data: dataSpec,
            params: paramsSpec
        ]

        // Add the stage to the output specification
        output_specification = output_specification + [(stageName): stageSpec]

        // Finally, store the hash for the next iteration
        prevHash = hash
    }

    // Process each organelle stage. If the stage is not enabled, we return nothing
    // Otherwise the hash is calculated using the parameters as well as the stage name,
    // to differentiate organelles and plastids.
    ORGANELLE_SPEC.each { stageName, config ->
        if(!config.enabled) return

        // Force organelle assembly mode "contigs"
        def mode = "contigs"

        // Organelle assembly always depends on the assembly stage
        prevHash = stageHashes["assembly"]

        // Extract data and params into submaps
        def dataSpec = spec.subMap(config.data)
        def paramsSpec = spec.subMap(config.params)

        // We add the stageName (mito or plastid) to the hash to ensure uniqueness
        // Concatenate the data and params values to hash and generate it
        // If we have a previous hash, we add this to the computed hash and recompute
        def hashContent = dataSpec.values().join("") + paramsSpec.values().join("") + stageName + mode
        def hash = (prevHash + hashContent).sha256()

        stageSpec = [
            hash: hash,
            prevHash: prevHash,
            data: dataSpec,
            params: paramsSpec + [organelle: stageName, mode: mode]
        ]

        output_specification = output_specification + [(stageName): stageSpec]
    }

    return output_specification
}

//
// Process a flat assembly specification for oatk
//
def stageOatkSpec(spec) {
    // Data required for oatk run
    def OATK_DATA = [
        "long_read_dataset",
        "long_read_platform"
    ]

    // Params required for oatk run
    def OATK_PARAMS = [
        "coverage",
        "oatk_kmer_size",
        "oatk_coverage_cutoff",
        "oatk_arguments",
        "oatk_mito_hmm",
        "oatk_plastid_hmm"
    ]

    // Extract data and params from spec
    def dataSpec = spec.subMap(OATK_DATA)
    def paramsSpec = spec.subMap(OATK_PARAMS)

    // Generate hash from configuration to detect changes
    def hash_input = dataSpec.values().join("") + paramsSpec.values().flatten().join("")
    def hash = hash_input.sha256()

    def stageSpec = [
        hash: hash,
        data: dataSpec,
        params: paramsSpec
    ]

    return spec.subMap(["id", "assembler"]) + [oatk: stageSpec]
}

//
// Process a flat assembly specification for mitohifi
//
def stageMitohifiSpec(spec) {
    // Data required for mitohifi run
    def MITOHIFI_DATA = [
        "long_read_dataset",
        "long_read_platform"
    ]

    // Params required for mitohifi run
    def MITOHIFI_PARAMS = [
        "mitohifi_reference_species",
        "mitohifi_mito_genetic_code",
        "mitohifi_arguments"
    ]

    // Extract data and params from spec
    def dataSpec = spec.subMap(MITOHIFI_DATA)
    def paramsSpec = spec.subMap(MITOHIFI_PARAMS)

    // Generate hash from configuration to detect changes
    def hash_input = dataSpec.values().join("") + paramsSpec.values().flatten().join("")
    def hash = hash_input.sha256()

    def stageSpec = [
        hash: hash,
        data: dataSpec,
        params: paramsSpec + [organelle: "mito", mode: "reads"]
    ]

    return spec.subMap(["id", "assembler"]) + [mitohifi: stageSpec]
}

//
// Add read file and FASTK dataset data to a spec .data list
//
def addData(spec, data_list, merqury_haptabs) {

    // Iterate over the key-value pairs in the spec and modify those
    // that are maps in place
    def specWithData = spec.collectEntries { key, value ->
        if (value instanceof Map) {
            // Pull out the dataset list
            def dataMap = value.data

            def allDataTypes = ["long_read", "ultralong", "hic", "polishing", "maternal", "paternal"]

            // Get the data types used
            def usedDataTypes = dataMap.keySet()
                .findAll { k -> k.endsWith("dataset") && dataMap[k] }
                .collect { type -> type - ~/_dataset/ }

            def updatedDataMap = allDataTypes.inject(dataMap) { data_map, dataType ->
                if(dataType in usedDataTypes) {
                    // Get dataset ID and platform
                    def dataID = data_map[(dataType + "_dataset")]
                    def dataPlatform = data_map[(dataType + "_platform")]

                    // Find dataset
                    def dataSet = data_list.find { data ->
                        data.id == dataID && data.platform == dataPlatform
                    }

                    data_map + [
                        (dataType + "_reads"): dataSet.reads,
                        (dataType + "_fk_hist"): dataSet.fk_hist,
                        (dataType + "_fk_ktab"): dataSet.fk_ktab,
                        (dataType + "_yak"): dataSet.yak,
                    ]
                } else {
                    data_map + [
                        (dataType + "_reads"): [],
                        (dataType + "_fk_hist"): [],
                        (dataType + "_fk_ktab"): [],
                        (dataType + "_yak"): [],
                    ]
                }
            }

            // Determine if this is a "trio" assembly
            def trioAssembly = ["maternal_dataset", "maternal_platform", "paternal_dataset", "paternal_platform"]
                .every { dataset -> dataMap[dataset] }

            // Add maternal and paternal haptabs to the dataset
            if(trioAssembly && merqury_haptabs) {
                def haptabs = merqury_haptabs.find { data ->
                    data.long_read_dataset == dataMap.long_read_dataset &&
                    data.platform == dataMap.long_read_platform &&
                    data.maternal_illumina_dataset == dataMap.illumina_maternal_dataset &&
                    data.paternal_illumina_dataset == dataMap.illumina_paternal_dataset
                }

                updatedDataMap = updatedDataMap + [
                    maternal_haptab: mat_haptab,
                    paternal_haptab: pat_haptab,
                ]
            } else {
                updatedDataMap = updatedDataMap + [
                    maternal_haptab: [],
                    paternal_haptab: [],
                ]
            }

            return [(key): value + [data: updatedDataMap]]
        } else {
            return [(key): value]
        }
    }

    return specWithData
}
