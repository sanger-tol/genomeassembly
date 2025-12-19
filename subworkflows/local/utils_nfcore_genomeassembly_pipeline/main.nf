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
include { READ_YAML                 } from '../../../modules/local/read_yaml'

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
    // Logic: read assembly specifications from the assembly specification sheet
    //        and validate various entries
    ch_assembly_specs = channel.fromList(
            samplesheetToList(asm_specs, "${projectDir}/assets/schema_assembly.json")
        )
        .map { list ->
            // Because a list with a single element is not destructured, extract spec explicitly
            def spec = list[0]

            // For validation, a list of all Hifiasm error correction arguments
            def hifiasm_error_correction_flags = [
                "-k", "-w", "-f", "-D", "-N", "-r", "-z", "--max-kocc", "--hg-size", "--min-hist-cnt"
            ]
            def mitohifi_disallowed_args = ["-c", "-r", "-f", "-g", "-t", "-o"]

            def invalid_hifiasm_ec_flags = validateArgString(spec.hifiasm_error_correction_arguments, hifiasm_error_correction_flags, false)
            def invalid_hifiasm_flags = validateArgString(spec.hifiasm_arguments, hifiasm_error_correction_flags, true)
            def invalid_oatk_args = validateArgString(spec.oatk_arguments, ["-k", "-c"], true)
            def invalid_mitohifi_contigs_args = validateArgString(spec.mitohifi_contigs_arguments, mitohifi_disallowed_args, true)
            def invalid_mitohifi_reads_args = validateArgString(spec.mitohifi_reads_arguments, mitohifi_disallowed_args, true)

            if(invalid_hifiasm_ec_flags) {
                error("Assembly specification error [${spec.id}]: hifiasm_error_correction_arguments contains invalid flag: ${invalid_hifiasm_ec_flags}")
            }

            if(invalid_hifiasm_flags) {
                error("Assembly specification error [${spec.id}]: hifiasm_arguments contains invalid flag: ${invalid_hifiasm_flags}")
            }

            if(invalid_oatk_args) {
                error("Assembly specification error [${spec.id}]: oatk_arguments contains invalid flag: ${invalid_oatk_args}")
            }

            if(invalid_mitohifi_contigs_args) {
                error("Assembly specification error [${spec.id}]: oatk_arguments contains invalid flag: ${invalid_mitohifi_contigs_args}")
            }

            if(invalid_mitohifi_reads_args) {
                error("Assembly specification error [${spec.id}]: oatk_arguments contains invalid flag: ${invalid_mitohifi_reads_args}")
            }

            if(spec.trio_assembly && (spec.maternal_illumina_dataset == spec.paternal_illumina_dataset)) {
                error("Assembly specification error [${spec.id}]: maternal and paternal Illumina datasets are the same")
            }

            if(spec.phased_assembly && (spec.long_read_dataset != spec.illumina_hic_dataset)) {
                log.warn("Assembly specification warning [${spec.id}]: Phased assembly is enabled, but the long read and Hi-C datasets are not the same. This might lead to incorrect results.")
            }

            return [spec]
        }

    //
    // Logic: read raw genomic data from the input sheet
    //
    ch_genomic_data = channel.fromList(
            samplesheetToList(genomic_data, "${projectDir}/assets/schema_genomic.json")
        )
        .collect().map { list -> [list] } 
        .combine(ch_assembly_specs)
        .multiMap { datasets, spec ->      
            // Validate and extract long read data
            def long_read_data = findDataset(spec.long_read_dataset, datasets)
            if (!long_read_data) {
                error("Assembly specification error [${spec.id}]: long read dataset '${spec.long_read_dataset}' not found")
            }
            
            def platform_key = spec.long_read_platform
            if (!hasData(long_read_data, platform_key)) {
                error("Assembly specification error [${spec.id}]: ${platform_key} data missing for dataset '${spec.long_read_dataset}'")
            }
            
            // Validate and extract Hi-C data (if scaffolding or phasing is enabled)
            def illumina_hic_data = null
            if ((spec.scaffold || spec.phased_assembly) && spec.illumina_hic_dataset) {
                illumina_hic_data = findDataset(spec.illumina_hic_dataset, datasets)
                if (!illumina_hic_data) {
                    error("Assembly specification error [${spec.id}]: Hi-C dataset '${spec.illumina_hic_dataset}' not found")
                }
                if (!hasData(illumina_hic_data, "illumina_hic")) {
                    error("Assembly specification error [${spec.id}]: illumina_hic data missing for dataset '${spec.illumina_hic_dataset}'")
                }
            }
            
            // Validate and extract 10X data (if specified)
            def illumina_10x_data = null
            if (spec.polishing && spec.illumina_10x_dataset) {
                illumina_10x_data = findDataset(spec.illumina_10x_dataset, datasets)
                if (!illumina_10x_data) {
                    error("Assembly specification error [${spec.id}]: 10X dataset '${spec.illumina_10x_dataset}' not found")
                }
                if (!hasData(illumina_10x_data, "illumina_10x")) {
                    error("Assembly specification error [${spec.id}]: illumina_10x data missing for dataset '${spec.illumina_10x_dataset}'")
                }
            }
            
            // Validate and extract trio data (if trio assembly is enabled)
            def maternal_illumina_data = null
            def paternal_illumina_data = null
            if (spec.trio_assembly) {
                if (!spec.maternal_illumina_dataset || !spec.paternal_illumina_dataset) {
                    error("Assembly specification error [${spec.id}]: trio_assembly enabled but maternal/paternal datasets not specified")
                }
                
                if (spec.maternal_illumina_dataset == spec.paternal_illumina_dataset) {
                    error("Assembly specification error [${spec.id}]: maternal and paternal Illumina datasets must be different")
                }
                
                maternal_illumina_data = findDataset(spec.maternal_illumina_dataset, datasets)
                paternal_illumina_data = findDataset(spec.paternal_illumina_dataset, datasets)
                
                if (!maternal_illumina_data) {
                    error("Assembly specification error [${spec.id}]: maternal dataset '${spec.maternal_illumina_dataset}' not found")
                }
                if (!paternal_illumina_data) {
                    error("Assembly specification error [${spec.id}]: paternal dataset '${spec.paternal_illumina_dataset}' not found")
                }
                
                if (!hasData(maternal_illumina_data, "illumina")) {
                    error("Assembly specification error [${spec.id}]: illumina data missing for maternal dataset '${spec.maternal_illumina_dataset}'")
                }
                if (!hasData(paternal_illumina_data, "illumina")) {
                    error("Assembly specification error [${spec.id}]: illumina data missing for paternal dataset '${spec.paternal_illumina_dataset}'")
                }
                
                // Check FastK kmer sizes match for trio
                def maternal_kmer = maternal_illumina_data.illumina?.fastk ? maternal_illumina_data.illumina.fastk[2] : null
                def paternal_kmer = paternal_illumina_data.illumina?.fastk ? paternal_illumina_data.illumina.fastk[2] : null
                def child_kmer = long_read_data[platform_key]?.fastk ? long_read_data[platform_key].fastk[2] : null
                
                if (maternal_kmer && paternal_kmer && child_kmer && 
                    [maternal_kmer, paternal_kmer, child_kmer].unique().size() > 1) {
                    error("Assembly specification error [${spec.id}]: FastK kmer sizes must match for trio (mat:${maternal_kmer}, pat:${paternal_kmer}, child:${child_kmer})")
                }
            }
            
            // Emit channels: [meta, reads, hist, ktab]
            long_reads: [
                [
                    id: spec.id,
                    kmer_size: long_read_data[platform_key].fastk ? long_read_data[platform_key].fastk[2] : null,
                    read_type: "long",
                    read_platform: spec.long_read_platform
                ],
                long_read_data[platform_key].reads ?: [],
                long_read_data[platform_key].fastk ? long_read_data[platform_key].fastk[0] : [],
                long_read_data[platform_key].fastk ? long_read_data[platform_key].fastk[1] : []
            ]
            
            illumina_hic: illumina_hic_data ? [
                [
                    id: spec.id,
                    read_type: "hic",
                ],
                illumina_hic_data.illumina_hic.reads ?: [],
                illumina_hic_data.illumina_hic.fastk ? illumina_hic_data.illumina_hic.fastk[0] : [],
                illumina_hic_data.illumina_hic.fastk ? illumina_hic_data.illumina_hic.fastk[1] : []
            ] : []
            
            illumina_10x: illumina_10x_data ? [
                [id: spec.id, read_type: "10x"],
                illumina_10x_data.illumina_10x.reads ?: [],
                illumina_10x_data.illumina_10x.fastk ? illumina_10x_data.illumina_10x.fastk[0] : [],
                illumina_10x_data.illumina_10x.fastk ? illumina_10x_data.illumina_10x.fastk[1] : []
            ] : []
            
            trio: (maternal_illumina_data && paternal_illumina_data) ? [
                [
                    [
                        id: spec.id,
                        kmer_size: maternal_illumina_data.illumina?.fastk ? maternal_illumina_data.illumina.fastk[2] : null,
                        trio: "mat"
                    ],
                    maternal_illumina_data.illumina.reads ?: [],
                    maternal_illumina_data.illumina?.fastk ? maternal_illumina_data.illumina.fastk[0] : [],
                    maternal_illumina_data.illumina?.fastk ? maternal_illumina_data.illumina.fastk[1] : []
                ],
                [
                    [
                        id: spec.id,
                        kmer_size: paternal_illumina_data.illumina?.fastk ? paternal_illumina_data.illumina.fastk[2] : null,
                        trio: "pat"
                    ],
                    paternal_illumina_data.illumina.reads ?: [],
                    paternal_illumina_data.illumina?.fastk ? paternal_illumina_data.illumina.fastk[0] : [],
                    paternal_illumina_data.illumina?.fastk ? paternal_illumina_data.illumina.fastk[1] : []
                ],
                [
                    [
                        id: spec.id,
                        kmer_size: long_read_data[platform_key].fastk ? long_read_data[platform_key].fastk[2] : null,
                        trio: "child"
                    ],
                    long_read_data[platform_key].reads ?: [],
                    long_read_data[platform_key].fastk ? long_read_data[platform_key].fastk[0] : [],
                    long_read_data[platform_key].fastk ? long_read_data[platform_key].fastk[1] : []
                ]
            ] : []
        }

    ch_genomic_data.illumina_10x.view()

    emit:
    specs         = ch_assembly_specs
    long_reads    = ch_genomic_data.long_reads.filter { _meta, reads, _hist, _ktab -> reads }.unique()
    hic_reads     = ch_genomic_data.illumina_hic.filter { _meta, reads -> reads }.unique()
    illumina_10x  = ch_genomic_data.illumina_10x.filter { _meta, reads -> reads }.unique()
    trio          = ch_genomic_data.trio.filter { _meta, mat, pat, child -> mat && pat && child }.unique()
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
def validateArgString(argString, flagList, rejectFlags) {
    def args = argString
        .tokenize()
        .collect { arg -> 
            if (arg =~ /^-[a-zA-Z]/) {
                // Single dash: extract just the flag letter
                arg.replaceAll(/^(-[a-zA-Z]).*/, '$1')
            } else if (arg =~ /^--[a-zA-Z]+/) {
                arg.replaceAll(/^(--[a-zA-Z]+).*/, '$1')
            } else {
                null
            }
        }
        .findAll { arg -> arg != null }
    
    // If rejectFlags, return the first arg found that is in the flag list
    // Otherwise, return the first arg found that is /not/ in the flag list
    if(rejectFlags) {
        return args.find { arg -> (arg in flagList) }
    } else {
        return args.find { arg -> !(arg in flagList) }
    }
}

// Check if dataset has data for the specified key
def hasData(dataMap, key) {
    return dataMap && dataMap[key] && dataMap[key].reads
}

// Find the dataset by ID
def findDataset(id, datasets) {
    if (!id) return null
    return datasets.find { data -> data.id == id }
}