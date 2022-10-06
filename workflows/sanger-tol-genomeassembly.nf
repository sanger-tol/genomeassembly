/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
//WorkflowSanger-tol-genomeassembly.initialise(params, log)

// Check input path parameters to see if they exist
def checkPathParamList = [ params.input ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = Channel.of(file(params.input)) } else { exit 1, 'Input samplesheet not specified!' }
if (params.groups) { groups = params.groups } else { groups = 100; }

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
include { ALIGN_SHORT   } from '../subworkflows/local/align_short'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/modules/custom/dumpsoftwareversions/main'
include { BWAMEM2_INDEX } from '../modules/nf-core/modules/bwamem2/index/main'

//
// MODULE: Installed from sanger-tol mirror nf-core/modules
//
include { LONGRANGER_MKREF } from '../modules/sanger-tol/nf-core-modules/longranger/mkref/main'
include { LONGRANGER_ALIGN } from '../modules/sanger-tol/nf-core-modules/longranger/align/main'
include { EXTRACT_SEQUENCES as EXTRACT_SEQUENCES_PRIMARY } from '../modules/local/extract_sequences'
include { EXTRACT_SEQUENCES as EXTRACT_SEQUENCES_HAPLOTIGS } from '../modules/local/extract_sequences'
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

workflow SANGER_TOL_GENOMEASSEMBLY {

    ch_versions = Channel.empty()

    //
    // SUBWORKFLOW: Read in yaml, validate and prepare for further steps
    //   
    PREPARE_INPUT(ch_input)
    ch_versions = ch_versions.mix(PREPARE_INPUT.out.versions)
    

    //
    // Polishing step 1: map reads to the reference
    //
    PREPARE_INPUT.out.assemblies.map{ meta, p, h, merged -> [meta, merged] }.set{ fasta_merged_ch }
    
    LONGRANGER_MKREF(fasta_merged_ch)
    ch_versions = ch_versions.mix(LONGRANGER_MKREF.out.versions)

    PREPARE_INPUT.out.illumina_10X.map{ meta, reads, kmers -> [reads]}
                    .set{ illumina_10X_ch }
    LONGRANGER_ALIGN( LONGRANGER_MKREF.out.folder, illumina_10X_ch )
    ch_versions = ch_versions.mix(LONGRANGER_ALIGN.out.versions)

    //
    // Polishing step 2: apply freebayes consensus based on longranger alignments
    //
    LONGRANGER_ALIGN.out.bam.join( LONGRANGER_ALIGN.out.bai ).set{ bam_ch }
    PREPARE_INPUT.out.assemblies.join( PREPARE_INPUT.out.indices )
                                .map{ meta, p, h, merged, p_i, h_i, merged_i -> [ merged, merged_i ] }
                                .set{ reference_ch }

    POLISHING(bam_ch, reference_ch, groups, LONGRANGER_ALIGN.out.csv.collect{it[1]} )    
    ch_versions = ch_versions.mix(POLISHING.out.versions)

    PREPARE_INPUT.out.indices.map{ meta, p_i, h_i, merged_i -> [p_i]}.set{primary_index_ch}
    EXTRACT_SEQUENCES_PRIMARY( POLISHING.out.fasta, primary_index_ch ) 
    PREPARE_INPUT.out.indices.map{ meta, p_i, h_i, merged_i -> [h_i]}.set{haplotigs_index_ch}
    EXTRACT_SEQUENCES_HAPLOTIGS( POLISHING.out.fasta, haplotigs_index_ch )

    PREPARE_INPUT.out.hic.view()
    PREPARE_INPUT.out.hic.map{ meta, crams, motif -> [meta, crams] }
                         .set{ crams_ch }

    BWAMEM2_INDEX ( EXTRACT_SEQUENCES_PRIMARY.out.subseq )
    ch_index = BWAMEM2_INDEX.out.index


    ALIGN_SHORT( crams_ch, ch_index, EXTRACT_SEQUENCES_PRIMARY.out.subseq.map{ meta, fasta -> [fasta]} )    
//    ALIGN_SHORT( Channel.of([[id:"ilEupCent1", datatype:"hic", read_group:"\'@RG\\tID:ilEupCent1\\tPL:ILLUMINA\\tSM:ilEupCent1\'"], "/lustre/scratch124/tol/projects/darwin/users/kk16/development/nextflow/sanger-tol-genomeassembly/work/de/7c601842368d43b1f06c765c67dbc3/ilEupCent1.markdup.bam"]), Channel.of([[id:"ilEupCent1"], "/lustre/scratch124/tol/projects/darwin/users/kk16/development/nextflow/sanger-tol-genomeassembly/work/ac/12febb43634cc9cc688d82916c49a2/bwamem2"]), Channel.of("/lustre/scratch124/tol/projects/darwin/users/kk16/development/nextflow/sanger-tol-genomeassembly/work/e5/2b9e6a9cb53eb48ff9c5a123cc5969/primary.fa"))    
    ALIGN_SHORT.out.bed.view() 
 
    
  
//    ch_bed = Channel.of([meta, bed])
//    ch_fa = Channel.of([fasta, fasta+".fai"])
//    SCAFFOLDING(ch_bed, ch_fa, true, '', '', 1000 )

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
