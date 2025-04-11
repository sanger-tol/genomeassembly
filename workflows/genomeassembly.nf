/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { KMERS         } from '../subworkflows/local/kmers'
include { POLISHING_10X } from '../subworkflows/local/polishing_10X'
include { PURGING       } from '../subworkflows/local/purging'
include { RAW_ASSEMBLY  } from '../subworkflows/local/raw_assembly'

//include { SCAFFOLDING                             } from '../subworkflows/local/scaffolding'
//include { ORGANELLES                              } from '../subworkflows/local/organelles'
//include { KEEP_SEQNAMES as KEEP_SEQNAMES_PRIMARY  } from '../modules/local/keep_seqnames'
//include { HIC_MAPPING                             } from '../subworkflows/local/hic_mapping'
//include { GENOME_STATISTICS } from '../subworkflows/local/genome_statistics'
//include { CAT_CAT as CAT_CAT_MITOHIFI_READS          } from "../modules/nf-core/cat/cat/main"
//include { SAMTOOLS_FAIDX as SAMTOOLS_FAIDX_PURGEDUPS } from '../modules/nf-core/samtools/faidx/main'
//include { SEQTK_SUBSEQ as SEQTK_SUBSEQ_HAPLOTIGS     } from '../modules/nf-core/seqtk/subseq/main'
//include { SEQTK_SUBSEQ as SEQTK_SUBSEQ_PRIMARY       } from '../modules/nf-core/seqtk/subseq/main'
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_genomeassembly_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow GENOMEASSEMBLY {

    take:
    long_reads
    hic_reads
    illumina_10x
    mat_reads
    pat_reads
    busco
    mito
    plastid

    main:

    ch_versions = Channel.empty()

    //
    // Subworkflow: generate Kmer databases and estimate
    //     coverage if not provided
    //
    KMERS(long_reads, mat_reads, pat_reads)
    ch_versions = ch_versions.mix(KMERS.out.versions)

    // FastK databases
    ch_fastk = KMERS.out.fastk

    // Get the long reads out with additional metadata from
    // the kmer-based analyses
    ch_long_reads = KMERS.out.long_reads
        | collect

    // Drop FastK databases for all downstream uses of Hi-C
    // CRAM files
    ch_hic_reads  = hic_reads
        | map { meta, cram, hist, ktab -> [meta, cram] }
        | collect

    ch_i10x_reads = illumina_10x
        | map { meta, reads, hist, ktab -> [meta, reads] }
        | collect

    ch_trio_yak_dbs   = KMERS.out.trio_yakdb

    //
    // Subworkflow: raw assembly of long reads using hifiasm
    //
    RAW_ASSEMBLY(
        ch_long_reads,
        ch_hic_reads,
        ch_trio_yak_dbs
    )
    ch_versions = ch_versions.mix(RAW_ASSEMBLY.out.versions)
    ch_assemblies = RAW_ASSEMBLY.out.assembly_fasta

    //
    // Subworkflow: Purge dups on specified assemblies
    //
    ch_assemblies_purge_status = ch_assemblies
        | branch { meta, assembly ->
            def purge_types = params.purging_assemblytypes.tokenize(",")
            purge: meta.assembly_type in purge_types
            no_purge: true
        }

    PURGING(
        ch_assemblies_purge_status.purge,
        ch_long_reads
    )
    ch_versions = ch_versions.mix(PURGING.out.versions)

    ch_assemblies = ch_assemblies
        | mix(PURGING.out.assemblies)

    ch_assemblies_to_polish = ch_assemblies
        | branch { meta, assembly ->
            def polish_types = params.polishing_assemblytypes.tokenize(",")
            polish: (params.enable_polishing && meta.assembly_type in polish_types && params.polishing_longranger_container_path)
            no_polish: true
        }

    POLISHING_10X(
        ch_assemblies_to_polish.polish,
        ch_i10x_reads
    )
    ch_versions = ch_versions.mix(POLISHING_10X.out.versions)

    ch_assemblies = ch_assemblies
        | mix(POLISHING_10X.out.assemblies)

//    //
//    // LOGIC: CREATE A CHANNEL FOR THE PATHS TO HIC DATA
//    //
//    hic.map{ meta, crams, motif, hic_aligner -> [meta, crams] }
//        .set{ crams_ch }
//
//    hic.map{ meta, crams, motif, hic_aligner -> [meta, hic_aligner] }
//        .set{ hic_aligner_ch }
//
//    //
//    // SUBWORKFLOW: MAP HIC DATA TO THE PRIMARY ASSEMBLY
//    //
//    HIC_MAPPING ( primary_contigs_ch,crams_ch,hic_aligner_ch, "")
//    ch_versions = ch_versions.mix(HIC_MAPPING.out.versions)
//
//    //
//    // SUBWORKFLOW: SCAFFOLD THE PRIMARY ASSEMBLY
//    //
//    SCAFFOLDING( HIC_MAPPING.out.bed, primary_contigs_ch, params.cool_bin, "")
//    ch_versions = ch_versions.mix(SCAFFOLDING.out.versions)
//
//    //
//    // LOGIC: CREATE A CHANNEL FOR THE FINAL ASSEMBLY REPRESENTED BY
//    //        THE SCAFFOLDS AND HAPLOTIGS
//    //
//    SCAFFOLDING.out.fasta.combine(haplotigs_ch)
//                        .map{meta_s, fasta_s, meta_h, fasta_h -> [ meta_h, fasta_s, fasta_h ]}
//                        .set{ stats_input_ch }
//
//    //
//    // SUBWORKFLOW: CALCULATE ASSEMBLY STATISTICS FOR THE FINAL ASSEMBLY
//    //
//    GENOME_STATISTICS_SCAFFOLDS( stats_input_ch,
//        busco,
//        GENOMESCOPE_MODEL.out.hist,
//        GENOMESCOPE_MODEL.out.ktab,
//        [],
//        [],
//        [],
//        [],
//        false
//    )
//
//    if ( params.hifiasm_hic_on ) {
//        //
//        // SUBWORKFLOW: MAP HIC DATA TO THE HAP1 CONTIGS
//        //
//        HIC_MAPPING_HAP1 ( RAW_ASSEMBLY.out.hap1_hic_contigs, crams_ch, hic_aligner_ch, 'hap1' )
//        ch_versions = ch_versions.mix(HIC_MAPPING_HAP1.out.versions)
//
//        //
//        // SUBWORKFLOW: SCAFFOLD HAP1
//        //
//        SCAFFOLDING_HAP1( HIC_MAPPING_HAP1.out.bed, RAW_ASSEMBLY.out.hap1_hic_contigs, params.cool_bin, 'hap1' )
//        ch_versions = ch_versions.mix(SCAFFOLDING_HAP1.out.versions)
//
//        //
//        // SUBWORKFLOW: MAP HIC DATA TO THE HAP2 CONTIGS
//        //
//        HIC_MAPPING_HAP2 ( RAW_ASSEMBLY.out.hap2_hic_contigs, crams_ch, hic_aligner_ch, 'hap2' )
//        ch_versions = ch_versions.mix(HIC_MAPPING_HAP2.out.versions)
//
//        //
//        // SUBWORKFLOW: SCAFFOLD HAP2
//        //
//        SCAFFOLDING_HAP2( HIC_MAPPING_HAP2.out.bed, RAW_ASSEMBLY.out.hap2_hic_contigs, params.cool_bin, 'hap2' )
//        ch_versions = ch_versions.mix(SCAFFOLDING_HAP2.out.versions)
//
//        //
//        // LOGIC: CREATE A CHANNEL FOR THE FULL HAP1/HAP2 ASSEMBLY
//        //
//        SCAFFOLDING_HAP1.out.fasta.combine(SCAFFOLDING_HAP2.out.fasta)
//                    .map{meta_s, fasta_s, meta_h, fasta_h -> [ [id:meta_h.id], fasta_s, fasta_h ]}
//                    .set{ stats_haps_input_ch }
//
//        //
//        // SUBWORKFLOW: CALCULATE ASSEMBLY STATISTICS FOR HAP1/HAP2 ASSEMBLY
//        //
//        GENOME_STATISTICS_SCAFFOLDS_HAPS( stats_haps_input_ch,
//            busco,
//            GENOMESCOPE_MODEL.out.hist,
//            GENOMESCOPE_MODEL.out.ktab,
//            [],
//            [],
//            [],
//            [],
//            true
//        )
//    }
//
//    if ( params.hifiasm_trio_on ) {
//        //
//        // SUBWORKFLOW: MAP HIC DATA TO THE PATERNAL CONTIGS
//        //
//        HIC_MAPPING_PAT ( RAW_ASSEMBLY.out.pat_contigs, crams_ch, hic_aligner_ch, 'pat' )
//        ch_versions = ch_versions.mix(HIC_MAPPING_PAT.out.versions)
//
//        //
//        // SUBWORKFLOW: SCAFFOLD PAT
//        //
//        SCAFFOLDING_PAT( HIC_MAPPING_PAT.out.bed, RAW_ASSEMBLY.out.pat_contigs, params.cool_bin, 'pat' )
//        ch_versions = ch_versions.mix(SCAFFOLDING_PAT.out.versions)
//
//        //
//        // SUBWORKFLOW: MAP HIC DATA TO THE MATERNAL CONTIGS
//        //
//        HIC_MAPPING_MAT ( RAW_ASSEMBLY.out.mat_contigs, crams_ch, hic_aligner_ch, 'mat' )
//        ch_versions = ch_versions.mix(HIC_MAPPING_MAT.out.versions)
//
//        //
//        // SUBWORKFLOW: SCAFFOLD MAT
//        //
//        SCAFFOLDING_MAT( HIC_MAPPING_MAT.out.bed, RAW_ASSEMBLY.out.mat_contigs, params.cool_bin, 'mat' )
//        ch_versions = ch_versions.mix(SCAFFOLDING_MAT.out.versions)
//
//        //
//        // LOGIC: CREATE A CHANNEL FOR THE FULL PAT/MAT ASSEMBLY
//        //
//        SCAFFOLDING_PAT.out.fasta.combine(SCAFFOLDING_MAT.out.fasta)
//            .map{meta_s, fasta_s, meta_h, fasta_h -> [ [id:meta_h.id], fasta_s, fasta_h ]}
//            .set{ stats_trio_input_ch }
//
//        //
//        // SUBWORKFLOW: CALCULATE ASSEMBLY STATISTICS FOR TRIO ASSEMBLY
//        //
//        GENOME_STATISTICS_SCAFFOLDS_TRIO( stats_trio_input_ch,
//            busco,
//            GENOMESCOPE_MODEL.out.hist,
//            GENOMESCOPE_MODEL.out.ktab,
//            fastk_pktab,
//            phapktab_ch,
//            fastk_mktab,
//            mhapktab_ch,
//            true
//        )
//
//    }


//
//    if ( params.organelles_on ) {
//        //
//        // LOGIC: CREATE CHANNEL FOR PRIMARY AND ALT CONTIGS
//        //
//        primary_contigs_ch.join(haplotigs_ch)
//            .map{ meta, pri, alt -> [meta, [pri, alt]]}
//            .set{ raw_pri_alt_ch }
//        //
//        // MODULE: MERGE PAW CONTIGS AND HAPLOTIGS INTO ONE FILE
//        //
//        CAT_CAT_RAW( raw_pri_alt_ch )
//
//        //
//        // LOGIC: DEFINE MERGED ASSEMBLY
//        //
//        merged_pri_alt_raw = CAT_CAT_RAW.out.file_out
//
//        //
//        // MODULE: MERGE INPUT FASTA FILES WITH PACBIO READS
//        //
//        CAT_CAT_MITOHIFI_READS(hifi_reads)
//        ch_versions = ch_versions.mix(CAT_CAT_MITOHIFI_READS.out.versions)
//
//        //
//        // SUBWORKFLOW: INDETIFY MITO IN THE RAW READS AND ASSEMBLY CONTIGS
//        //
//        ORGANELLES(CAT_CAT_MITOHIFI_READS.out.file_out, merged_pri_alt_raw,
//            mito, plastid)
//        ch_versions = ch_versions.mix(ORGANELLES.out.versions)
//    }
//

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name:  'genomeassembly_software_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }


    emit:
    versions       = ch_versions                 // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
