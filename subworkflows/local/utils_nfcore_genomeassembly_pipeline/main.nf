//
// Subworkflow with functionality specific to the sanger-tol/genomeassembly pipeline
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { UTILS_NFSCHEMA_PLUGIN     } from '../../nf-core/utils_nfschema_plugin'
include { paramsSummaryMap          } from 'plugin/nf-schema'
include { samplesheetToList         } from 'plugin/nf-schema'
include { paramsHelp                } from 'plugin/nf-schema'
include { completionEmail           } from '../../nf-core/utils_nfcore_pipeline'
include { completionSummary         } from '../../nf-core/utils_nfcore_pipeline'
include { imNotification            } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NFCORE_PIPELINE     } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NEXTFLOW_PIPELINE   } from '../../nf-core/utils_nextflow_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW TO INITIALISE PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_INITIALISATION {

    take:
    version           // boolean: Display version and exit
    validate_params   // boolean: Boolean whether to validate parameters against the schema at runtime
    monochrome_logs   // boolean: Do not use coloured log outputs
    nextflow_cli_args //   array: List of positional nextflow CLI args
    outdir            //  string: The output directory where the results will be saved
    help              // boolean: Display help message and exit
    help_full         // boolean: Show the full help message
    show_hidden       // boolean: Show hidden parameters in the help message
    genomic_data      //  string: Path to genome data samplesheet
    asm_specs         //  string: Path to assembly specification samplesheet

    main:
    //
    // Print version and exit if required and dump pipeline parameters to JSON file
    //
    UTILS_NEXTFLOW_PIPELINE (
        version,
        true,
        outdir,
        workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1
    )

    //
    // Validate parameters and generate parameter summary to stdout
    //
    before_text = """
-\033[2m----------------------------------------------------\033[0m-
\033[0;34m   _____                               \033[0;32m _______   \033[0;31m _\033[0m
\033[0;34m  / ____|                              \033[0;32m|__   __|  \033[0;31m| |\033[0m
\033[0;34m | (___   __ _ _ __   __ _  ___ _ __ \033[0m ___ \033[0;32m| |\033[0;33m ___ \033[0;31m| |\033[0m
\033[0;34m  \\___ \\ / _` | '_ \\ / _` |/ _ \\ '__|\033[0m|___|\033[0;32m| |\033[0;33m/ _ \\\033[0;31m| |\033[0m
\033[0;34m  ____) | (_| | | | | (_| |  __/ |        \033[0;32m| |\033[0;33m (_) \033[0;31m| |____\033[0m
\033[0;34m |_____/ \\__,_|_| |_|\\__, |\\___|_|        \033[0;32m|_|\033[0;33m\\___/\033[0;31m|______|\033[0m
\033[0;34m                      __/ |\033[0m
\033[0;34m                     |___/\033[0m
\033[0;35m  ${workflow.manifest.name} ${workflow.manifest.version}\033[0m
-\033[2m----------------------------------------------------\033[0m-
        """
    after_text = """${workflow.manifest.doi ? "\n* The pipeline\n" : ""}${workflow.manifest.doi.tokenize(",").collect { doi -> "    https://doi.org/${doi.trim().replace('https://doi.org/', '')}" }.join("\n")}${workflow.manifest.doi ? "\n" : ""}
* The nf-core framework
    https://doi.org/10.1038/s41587-020-0439-x

* Software dependencies
    https://github.com/nf-core/genomeassembly/blob/main/CITATIONS.md
"""
    command = "nextflow run ${workflow.manifest.name} -profile <docker/singularity/.../institute> --input samplesheet.csv --outdir <OUTDIR>"

    UTILS_NFSCHEMA_PLUGIN (
        workflow,
        validate_params,
        null,
        help,
        help_full,
        show_hidden,
        before_text,
        after_text,
        command
    )

    //
    // Check config provided to the pipeline
    //
    UTILS_NFCORE_PIPELINE (
        nextflow_cli_args
    )

    //
    // Logic: read raw genomic data from the input sheet and validate the dataset inputs.
    //
    ch_raw_genomic_data = channel.fromList(
            samplesheetToList(genomic_data, "${projectDir}/assets/schema_genomic.json")
        )
        .map { meta, reads, fastk ->
            // Check that all files have the same extension
            def extensions = reads.collect { file ->
                if(file.getExtension() == "gz") {
                    return file.getName().tokenize(".").takeLast(2).join('.')
                } else {
                    return file.getExtension()
                }
            }.unique()

            if(extensions.size() > 1) {
                error("Dataset validation error [${meta.id}:${meta.platform}]: Not all files are of the same type!")
            }

            // Validate that the correct extension is present for each platform
            if(meta.platform == "pacbio_hifi") {
                validateFileExtension(meta, extensions[0], /^fn?(asta)?(\.gz)?$/)
            } else if (meta.platform == "oxford_nanopore" || meta.platform == "illumina_10x") {
                validateFileExtension(meta, extensions[0], /^f(ast)?q(\.gz)?$/)
            } else {
                validateFileExtension(meta, extensions[0], /^cram$/)
            }

            // Construct a hybrid dataset name using the sample id and platform
            def datasetName = "${meta.id}.${meta.platform}"

            return [ meta + [id: datasetName], reads, fastk ]
        }

    //
    // Logic: read assembly specifications from the assembly specification sheet, and
    // validate that the requested data entries exist and that the optional tuning arguments
    // are valid.
    //
    ch_assembly_specs = channel.fromList(
            samplesheetToList(asm_specs, "${projectDir}/assets/schema_assembly.json")
        )
        // We need to wrap these up inside a list to have a single entry after combining
        .combine(ch_raw_genomic_data.collect().map { list -> [list] } )
        .map { spec, data ->
            // Add data type suffix to each dataset, e.g. "test" -> "test.illumina_hic"
            // This allows joining with the data downstream based on a single key
            spec = addDatasetSuffices(spec)

            // Check that each requested dataset exists in the input data schema
            checkDataExists(spec, data)

            // Validate any provided hifiasm args
            validateHifiasmArgs(spec.id, spec.hifiasm_error_correction_arguments, spec.hifiasm_arguments)

            // Validate trio assembly inputs
            if(spec.trio_assembly && (spec.maternal_illumina_dataset == spec.paternal_illumina_dataset)) {
                error("Assembly specification error [${spec.id}]: maternal and paternal Illumina datasets are the same")
            }

            // Validate phased assembly inputs
            if(spec.phased_assembly && (spec.long_read_dataset != spec.illumina_hic_dataset)) {
                log.warn("Assembly specification warning [${spec.id}]: Phased assembly is enabled, but the long read and Hi-C datasets are not the same. This might lead to incorrect results.")
            }

            return spec
        }

    //
    // Logic: After validation, create a channel containing the names of all the used
    // datasets to filter the whole dataset channel.
    //
    ch_used_datasets = ch_assembly_specs
        .map { spec ->
            return [
                spec.long_read_dataset,
                spec.illumina_hic_dataset,
                spec.illumina_10x_dataset,
                spec.maternal_illumina_dataset,
                spec.paternal_illumina_dataset
            ].findAll()
        }.transpose().unique().collect()

    //
    // Logic: Filter the input genomic datasets to only include those used in an assembly
    // specification
    //
    ch_genomic_data = ch_raw_genomic_data
        .combine(ch_used_datasets)
        .filter { meta, reads, fastk, data_list ->
            meta.id in data_list
        }
        .map { meta, reads, fastk, data_list -> [ meta, reads, fastk ] }
        .branch { meta, reads, fastk ->
            long_reads: meta.platform in ["pacbio_hifi", "oxford_nanopore"]
                return [ meta, reads, fastk ]
            illumina_hic: meta.platform == "illumina_hic"
                return [ meta, reads ]
            illumina_10x: meta.platform == "illumina_10x"
                return [ meta, reads ]
            illumina: meta.platform == "illumina"
                return [ meta, reads, fastk ]
        }

    emit:
    specs         = ch_assembly_specs
    long_reads    = ch_out_genomic_data.long_reads
    hic_reads     = ch_out_genomic_data.illumina_hic
    illumina_10x  = ch_out_genomic_data.illumina_10x
    illumina      = ch_out_genomic_data.illumina
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW FOR PIPELINE COMPLETION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_COMPLETION {

    take:
    email           //  string: email address
    email_on_fail   //  string: email address sent on pipeline failure
    plaintext_email // boolean: Send plain-text email instead of HTML
    outdir          //    path: Path to output directory where results will be published
    monochrome_logs // boolean: Disable ANSI colour codes in log output
    hook_url        //  string: hook URL for notifications

    main:
    summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")

    //
    // Completion email and summary
    //
    workflow.onComplete {
        if (email || email_on_fail) {
            completionEmail(
                summary_params,
                email,
                email_on_fail,
                plaintext_email,
                outdir,
                monochrome_logs,
                []
            )
        }

        completionSummary(monochrome_logs)
        if (hook_url) {
            imNotification(summary_params, hook_url)
        }
    }

    workflow.onError {
        log.error "Pipeline failed. Please refer to troubleshooting docs: https://nf-co.re/docs/usage/troubleshooting"
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// Validate channels from input samplesheet
//
def validateInputSamplesheet(input) {
    def (metas, fastqs) = input[1..2]

    // Check that multiple runs of the same sample are of the same datatype i.e. single-end / paired-end
    def endedness_ok = metas.collect{ meta -> meta.single_end }.unique().size == 1
    if (!endedness_ok) {
        error("Please check input samplesheet -> Multiple runs of a sample must be of the same datatype i.e. single-end or paired-end: ${metas[0].id}")
    }

    return [ metas[0], fastqs ]
}
//
// Generate methods description for MultiQC
//
def toolCitationText() {
    // Can use ternary operators to dynamically construct based conditions, e.g. params["run_xyz"] ? "Tool (Foo et al. 2023)" : "",
    // Uncomment function in methodsDescriptionText to render in MultiQC report
    def citation_text = [
            "Tools used in the workflow included:",
            "."
        ].join(' ').trim()

    return citation_text
}

def toolBibliographyText() {
    // Can use ternary operators to dynamically construct based conditions, e.g. params["run_xyz"] ? "<li>Author (2023) Pub name, Journal, DOI</li>" : "",
    // Uncomment function in methodsDescriptionText to render in MultiQC report
    def reference_text = [
        ].join(' ').trim()

    return reference_text
}

def methodsDescriptionText(mqc_methods_yaml) {
    // Convert  to a named map so can be used as with familiar NXF ${workflow} variable syntax in the MultiQC YML file
    def meta = [:]
    meta.workflow = workflow.toMap()
    meta["manifest_map"] = workflow.manifest.toMap()

    // Pipeline DOI
    if (meta.manifest_map.doi) {
        // Using a loop to handle multiple DOIs
        // Removing `https://doi.org/` to handle pipelines using DOIs vs DOI resolvers
        // Removing ` ` since the manifest.doi is a string and not a proper list
        def temp_doi_ref = ""
        def manifest_doi = meta.manifest_map.doi.tokenize(",")
        manifest_doi.each { doi_ref ->
            temp_doi_ref += "(doi: <a href=\'https://doi.org/${doi_ref.replace("https://doi.org/", "").replace(" ", "")}\'>${doi_ref.replace("https://doi.org/", "").replace(" ", "")}</a>), "
        }
        meta["doi_text"] = temp_doi_ref.substring(0, temp_doi_ref.length() - 2)
    } else meta["doi_text"] = ""
    meta["nodoi_text"] = meta.manifest_map.doi ? "" : "<li>If available, make sure to update the text to include the Zenodo DOI of version of the pipeline used. </li>"

    // Tool references
    meta["tool_citations"] = ""
    meta["tool_bibliography"] = ""

    // meta["tool_citations"] = toolCitationText().replaceAll(", \\.", ".").replaceAll("\\. \\.", ".").replaceAll(", \\.", ".")
    // meta["tool_bibliography"] = toolBibliographyText()


    def methods_text = mqc_methods_yaml.text

    def engine =  new groovy.text.SimpleTemplateEngine()
    def description_html = engine.createTemplate(methods_text).make(meta)

    return description_html.toString()
}

//
// sanger-tol/genomeassembly specific functions
//

// Validate a file extension against a regular expression and error if they do not match.
def validateFileExtension(meta, extension, regex) {
    if(!(extension =~ regex)) {
        error("Dataset validation error [${meta.id}:${meta.platform}]: File extension ${extension} does not match the expected input ${regex}.")
    }
}

// Validate a string of command line arguments against a list of flags.
//
// if rejectFlags is true, check if any flags in the string is in the flag list and return the first one.
// Otherwise, check if any flags in the string is /not/ in the flag list.
def validateArgString(argString, flagList, rejectFlags) {
    def args = argString
        .tokenize()
        .collect { arg ->
            def match = arg =~ /^(-{1,2}[a-zA-Z]+)(?=[^a-zA-Z]|$)/
            match ? match[0][0] : null
        }
        .findAll { it }

    return args.find { arg -> rejectFlags ? (arg in flagList) : !(arg in flagList) }
}

// Validate provided Hifiasm arguments to ensure that EC arguments are only specified in the
// hifiasm_error_correction_arguments assembly param
def validateHifiasmArgs(id, hifiasm_ec_args, hifiasm_args) {
    def valid_hifiasm_ec_flags: ["-k", "-w", "-f", "-D", "-N", "-r", "-z", "--max-kocc", "--hg-size", "--min-hist-cnt"]

    // Validate the hifiasm error-correction and normal arguments against the valid EC flag list
    [
        [hifiasm_ec_args, "hifiasm_error_correction_arguments", false],
        [hifiasm_args, "hifiasm_arguments", true]
    ].each { argString, argName, rejectFlags ->
        def invalid = validateArgString(argString, valid_hifiasm_ec_flags, rejectFlags)
        if(invalid) {
            error("Assembly specification error [${id}]: ${argName} contains invalid flag: ${invalid}")
        }
    }
}

// Add the dataset platform to the end of the dataset name to allow easy matching on a single key
//
// e.g. [illumina_hic_dataset: test1] -> [illumina_hic_dataset: test1.illumina_hic]
def addDatasetSuffices(spec) {
    // Map dataset names to platforms
    def datasetPlatforms = [
        long_read_dataset: spec.long_read_platform,
        illumina_hic_dataset: 'illumina_hic',
        illumina_10x_dataset: 'illumina_10x',
        maternal_illumina_dataset: 'illumina',
        paternal_illumina_dataset: 'illumina'
    ]

    // Rename datasets with platform suffixes
    return spec + datasetPlatforms.collectEntries { key, platform ->
        def value = spec[key]
        [(key): value ? "${value}.${platform}" : null]
    }
}

// For each dataset in an assembly specification, check that it exists in the data
// specification.
def checkDataExists(spec, data) {
    def datasetNames = [
        'long_read_dataset',
        'illumina_hic_dataset',
        'illumina_10x_dataset',
        'maternal_illumina_dataset',
        'paternal_illumina_dataset'
    ]

    datasetNames.each { datasetName ->
        def datasetID = spec[datasetName]
        if(datasetId && !data.find { it.id == datasetId }) {
            error("Assembly specification error [${spec.id}]: ${datasetName} '${datasetId}' does not exist!")
        }
    }
}
