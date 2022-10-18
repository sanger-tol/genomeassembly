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

if (params.motif) { motif = params.motif } else { motif = ''; }

if (params.resolutions) { resolutions = params.resolutions } else { resolutions = ''; }

if (params.cool_bin) { cool_bin = params.cool_bin } else { cool_bin = 1000; }

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
include { SCAFFOLDING   } from '../subworkflows/local/scaffolding'
include { GENOME_STATISTICS as GENOME_STATISTICS_POLISHED  } from '../subworkflows/local/assembly_stats'
include { GENOME_STATISTICS as GENOME_STATISTICS_SCAFFOLDED } from '../subworkflows/local/assembly_stats'

include { EXTRACT_SEQUENCES as EXTRACT_SEQUENCES_PRIMARY } from '../modules/local/extract_sequences'
include { EXTRACT_SEQUENCES as EXTRACT_SEQUENCES_HAPLOTIGS } from '../modules/local/extract_sequences'
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
include { SAMTOOLS_FAIDX } from '../modules/nf-core/modules/samtools/faidx/main'

//
// MODULE: Installed from sanger-tol mirror nf-core/modules
//
include { LONGRANGER_MKREF } from '../modules/sanger-tol/nf-core-modules/longranger/mkref/main'
include { LONGRANGER_ALIGN } from '../modules/sanger-tol/nf-core-modules/longranger/align/main'

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
    POLISHING.out.fasta.map{ meta, fasta -> [ [ id : meta.id, prefix: 'primary'], fasta ]}
                        .set{ primary_ch }
    EXTRACT_SEQUENCES_PRIMARY( primary_ch, primary_index_ch ) 
    ch_versions = ch_versions.mix(EXTRACT_SEQUENCES_PRIMARY.out.versions)

    PREPARE_INPUT.out.indices.map{ meta, p_i, h_i, merged_i -> [h_i]}.set{haplotigs_index_ch}
    POLISHING.out.fasta.map{ meta, fasta -> [ [ id : meta.id, prefix: 'haplotigs'], fasta ]}
                        .set{ haplotigs_ch }
    EXTRACT_SEQUENCES_HAPLOTIGS( haplotigs_ch, haplotigs_index_ch )
    ch_versions = ch_versions.mix(EXTRACT_SEQUENCES_HAPLOTIGS.out.versions)


    GENOME_STATISTICS_POLISHED( EXTRACT_SEQUENCES_PRIMARY.out.subseq.join(EXTRACT_SEQUENCES_HAPLOTIGS.out.subseq), 
                       PREPARE_INPUT.out.busco.map{ meta, buscodb, lineage -> buscodb},
                       PREPARE_INPUT.out.hifi.map{ meta, reads, kmerdb -> [meta, kmerdb]} )
    ch_versions = ch_versions.mix(GENOME_STATISTICS_POLISHED.out.versions)

    PREPARE_INPUT.out.hic.map{ meta, crams, motif -> [meta, crams] }
                         .set{ crams_ch }

    BWAMEM2_INDEX ( EXTRACT_SEQUENCES_PRIMARY.out.subseq )
    ch_versions = ch_versions.mix(BWAMEM2_INDEX.out.versions)
    ch_index = BWAMEM2_INDEX.out.index

    ALIGN_SHORT( crams_ch, ch_index, EXTRACT_SEQUENCES_PRIMARY.out.subseq.map{ meta, fasta -> [fasta]} )    
    ch_versions = ch_versions.mix(ALIGN_SHORT.out.versions)

    SAMTOOLS_FAIDX(EXTRACT_SEQUENCES_PRIMARY.out.subseq)
    ch_versions = ch_versions.mix(SAMTOOLS_FAIDX.out.versions)

    SAMTOOLS_FAIDX.out.fai.join(EXTRACT_SEQUENCES_PRIMARY.out.subseq)
                    .map{ meta, fai, fasta -> [fasta, fai] }
                    .set{ scaf_ref_ch }  
    
    SCAFFOLDING(ALIGN_SHORT.out.bed, scaf_ref_ch, true, motif, resolutions, cool_bin )
    ch_versions = ch_versions.mix(SCAFFOLDING.out.versions)
    
    GENOME_STATISTICS_SCAFFOLDED( SCAFFOLDING.out.fasta.join( EXTRACT_SEQUENCES_HAPLOTIGS.out.subseq ),
                                    PREPARE_INPUT.out.busco.map{ meta, buscodb, lineage -> buscodb},
                                    PREPARE_INPUT.out.hifi.map{ meta, reads, kmerdb -> [meta, kmerdb]})
    ch_versions = ch_versions.mix(GENOME_STATISTICS_SCAFFOLDED.out.versions)

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
