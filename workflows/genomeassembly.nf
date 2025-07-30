/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { HIC_MAPPING       } from '../subworkflows/sanger-tol/hic_mapping'
include { HIC_MAPPING_STATS } from '../subworkflows/local/hic_mapping_stats'
include { KMERS             } from '../subworkflows/local/kmers'
include { POLISHING_10X     } from '../subworkflows/local/polishing_10X'
include { PURGING           } from '../subworkflows/local/purging'
include { RAW_ASSEMBLY      } from '../subworkflows/local/raw_assembly'

//include { SCAFFOLDING                             } from '../subworkflows/local/scaffolding'
//include { ORGANELLES                              } from '../subworkflows/local/organelles'
//include { KEEP_SEQNAMES as KEEP_SEQNAMES_PRIMARY  } from '../modules/local/keep_seqnames'
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
        | map { meta, cram, _hist, _ktab -> [meta, cram] }
        | collect

    ch_i10x_reads = illumina_10x
        | map { meta, reads, _hist, _ktab -> [meta, reads] }
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

    //
    // Logic: Identify what steps to run on each assembly, using the params
    //
    ch_assemblies = RAW_ASSEMBLY.out.assembly_fasta
        | map { meta, assembly ->
            def purge_types  = params.purging_assemblytypes.tokenize(",")
            def polish_types = params.polishing_assemblytypes.tokenize(",")

            def purge        = meta.assembly_type in purge_types
            def polish       = (meta.assembly_type in polish_types && params.enable_polishing && params.polishing_longranger_container_path)
            def meta_new     = meta + [purge: purge, polish: polish]

            [meta_new, assembly]
        }

    //
    // Subworkflow: Purge dups on specified assemblies
    //
    ch_assemblies_to_purge = ch_assemblies
        | filter { meta, _assembly ->
            meta.purge == true
        }

    PURGING(
        ch_assemblies_to_purge,
        ch_long_reads
    )
    ch_versions = ch_versions.mix(PURGING.out.versions)

    //
    // Logic: Add metadata to purged assemblies
    //
    ch_purged_assemblies = PURGING.out.assemblies
        | map { meta, assembly ->
            def meta_new = meta + [assembly_stage: "purged"]
            [meta_new, assembly]
        }

    ch_assemblies = ch_assemblies
        | mix(ch_purged_assemblies)

    //
    // Logic: Filter the input assemblies so that if purging is enabled for an assembly type,
    //        only purged assemblies are polished.
    //
    ch_assemblies_to_polish = ch_assemblies
        | filter { meta, _assembly ->
            // If we have purged this type of assembly, remove raw stages
            def purging_filter  = (meta.purge == true) ? (meta.assembly_stage != "raw") : meta.assembly_stage == "raw"
            def polish_filter   = (meta.polish == true)
            def polishing_enabled = (params.enable_polishing && params.polishing_longranger_container_path)

            (polishing_enabled && polish_filter && purging_filter)
        }

    POLISHING_10X(
        ch_assemblies_to_polish,
        ch_i10x_reads
    )
    ch_versions = ch_versions.mix(POLISHING_10X.out.versions)

    //
    // Logic: Add metadata to polished assemblies
    //
    ch_polished_assemblies = POLISHING_10X.out.assemblies
        | map { meta, fasta ->
            [meta + [assembly_stage: "polished"], fasta]
        }

    ch_assemblies = ch_assemblies
        | mix(ch_polished_assemblies)

    //
    // Logic: Filter the assemblies to keep those for hi-c mapping
    //
    ch_assemblies_for_hic_mapping = ch_assemblies
        | filter { meta, _assembly ->
            if(meta.assembly_stage == "polished") { return true }
            else if(meta.assembly_stage == "purged" && !meta.polish) { return true }
            else if(meta.assembly_stage == "raw" && !meta.purge && !meta.polish) { return true }
            else { return false }
        }

    //
    // Subworkflow: Map Hi-C data to each assembly
    //
    HIC_MAPPING(
        ch_assemblies_for_hic_mapping,
        ch_hic_reads,
        params.hic_aligner,
        params.hic_mapping_cram_chunk_size,
        true // mark duplicates
    )
    ch_versions = ch_versions.mix(HIC_MAPPING.out.versions)

    //
    // Subworkflow: Calculate stats for Hi-C mapping
    //
    HIC_MAPPING_STATS(
        HIC_MAPPING.out.bam,
        ch_assemblies_for_hic_mapping
    )
    ch_versions = ch_versions.mix(HIC_MAPPING_STATS.out.versions)

    //
    // Subworkflow: scaffold assemblies
    //
    SCAFFOLDING(
        ch_assemblies_for_hic_mapping,
        HIC_MAPPING.out.bam,
        params.cool_bin
    )
    ch_versions = ch_versions.mix(SCAFFOLDING.out.versions)
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
