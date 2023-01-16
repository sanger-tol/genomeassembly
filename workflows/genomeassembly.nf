/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Check input path parameters to see if they exist
def checkPathParamList = [ params.input ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = Channel.of(file(params.input)) } else { exit 1, 'Input samplesheet not specified!' }
if (params.groups) { groups = params.groups } else { groups = 100; }

if (params.cool_bin) { cool_bin = params.cool_bin } else { cool_bin = 1000; }

if (params.polishing_on) { polishing_on = params.polishing_on } else { polishing_on = false; }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { PREPARE_INPUT } from '../subworkflows/local/prepare_input'
include { POLISHING     } from '../subworkflows/local/polishing'
include { SCAFFOLDING   } from '../subworkflows/local/scaffolding'
include { KEEP_SEQNAMES } from '../modules/local/keep_seqnames'
include { ALIGN_SHORT   } from '../subworkflows/local/align_short'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { SAMTOOLS_FAIDX   } from '../modules/nf-core/samtools/faidx/main.nf'

//
// MODULE: Installed directly from nf-core/modules
//
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/custom/dumpsoftwareversions/main'
include { SEQTK_SUBSEQ } from '../modules/nf-core/seqtk/subseq/main'
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

workflow GENOMEASSEMBLY {

    ch_versions = Channel.empty()

    //
    // SUBWORKFLOW: Read in yaml, validate and prepare for further steps
    //   
    PREPARE_INPUT(ch_input)
    ch_versions = ch_versions.mix(PREPARE_INPUT.out.versions)

    PREPARE_INPUT.out.assemblies.map{ meta, p, h, merged -> [meta, p] }.set{ primary_contigs_ch } 
    PREPARE_INPUT.out.assemblies.map{ meta, p, h, merged -> [meta, h] }.set{ haplotigs_ch } 

    if ( polishing_on ) {
        PREPARE_INPUT.out.illumina_10X.map{ meta, reads, kmers -> [reads]}
                        .set{ illumina_10X_ch }
        PREPARE_INPUT.out.assemblies.join(PREPARE_INPUT.out.indices)
                                    .map{ meta, p, h, merged, p_i, h_i, merged_i -> [ meta, merged, merged_i ] }
                                    .set{ reference_ch }
        POLISHING(reference_ch, illumina_10X_ch, groups)    
        ch_versions = ch_versions.mix(POLISHING.out.versions)

        // Separate the primary and alternative contigs again after polishing
        PREPARE_INPUT.out.indices.map{ meta, p_i, h_i, merged_i -> [meta, p_i]}.set{primary_index_ch}
        KEEP_SEQNAMES(primary_index_ch)
        ch_versions = ch_versions.mix(KEEP_SEQNAMES.out.versions)
        POLISHING.out.fasta.map{ meta, f -> f }
                           .set{ polished_fasta }
        SEQTK_SUBSEQ(polished_fasta, KEEP_SEQNAMES.out.seqlist)
        POLISHING.out.fasta.map{ meta, f -> meta }
                            .combine(SEQTK_SUBSEQ.out.sequences)
                            .set{ primary_contigs_ch }
    }
    
    PREPARE_INPUT.out.hic.map{ meta, crams, motif -> [meta, crams] }
                         .set{ crams_ch }

    // Map HiC data to the primary assembly
    primary_contigs_ch.map{ meta, fasta -> [ fasta ] }
                      .set{ hic_ref_ch }
    ALIGN_SHORT( crams_ch, hic_ref_ch )    
    ch_versions = ch_versions.mix(ALIGN_SHORT.out.versions)

    SCAFFOLDING( ALIGN_SHORT.out.bed, primary_contigs_ch, cool_bin )
    ch_versions = ch_versions.mix(SCAFFOLDING.out.versions)

    SCAFFOLDING.out.fasta.combine(haplotigs_ch)
                        .map{meta_s, fasta_s, meta_h, fasta_h -> [ meta_h, fasta_s, fasta_h ]}
                        .set{ stats_input_ch }

    //
    // MODULE: Collate versions.yml file
    //
    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.summary(workflow, params, log)
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
