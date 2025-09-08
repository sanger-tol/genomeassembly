/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { GENOME_STATISTICS } from '../subworkflows/sanger-tol/genome_statistics'
include { KMERS             } from '../subworkflows/local/kmers'
include { POLISHING_10X     } from '../subworkflows/local/polishing_10X'
include { PURGING           } from '../subworkflows/local/purging'
include { RAW_ASSEMBLY      } from '../subworkflows/local/raw_assembly'
include { SCAFFOLDING       } from '../subworkflows/local/scaffolding'

include { CAT_CAT as CONCATENATE_ASSEMBLIES     } from '../modules/nf-core/cat/cat'
include { SEQKIT_GREP as SEQKIT_GREP_SPLIT_HAPS } from '../modules/nf-core/seqkit/grep/main'


//include { ORGANELLES                              } from '../subworkflows/local/organelles'
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
    busco_lineage
    busco_lineage_directory
    mito_hmm
    plastid_hmm

    main:
    ch_versions = Channel.empty()

    //
    // Subworkflow: generate Kmer databases and estimate
    //     coverage if not provided
    //
    KMERS(long_reads, mat_reads, pat_reads)
    ch_versions = ch_versions.mix(KMERS.out.versions)

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

    //
    // Subworkflow: raw assembly of long reads using hifiasm
    //
    RAW_ASSEMBLY(
        ch_long_reads,
        ch_hic_reads,
        KMERS.out.trio_yakdb
    )
    ch_versions = ch_versions.mix(RAW_ASSEMBLY.out.versions)

    //
    // Logic: Identify what steps to run on each assembly, using the params
    //
    ch_assemblies_raw = RAW_ASSEMBLY.out.assembly_fasta
        | map { meta, hap1, hap2 ->
            def purge_types  = params.purging_assemblytypes.tokenize(",")
            def polish_types = params.polishing_assemblytypes.tokenize(",")

            def purge        = meta.assembly_type in purge_types
            def polish       = (meta.assembly_type in polish_types && params.enable_polishing && params.polishing_longranger_container_path)
            def meta_new     = meta + [purge: purge, polish: polish]

            [meta_new, hap1, hap2]
        }

    //
    // Subworkflow: Purge dups on specified assemblies
    //
    ch_assemblies_to_purge = ch_assemblies_raw
        | branch { meta, _hap1, _hap2 ->
            purge: meta.purge == true
            no_purge: true
        }

    PURGING(
        ch_assemblies_to_purge.purge,
        ch_long_reads
    )
    ch_versions = ch_versions.mix(PURGING.out.versions)

    //
    // Logic: Add metadata to purged assemblies
    //
    ch_purged_assemblies = PURGING.out.assemblies
        | map { meta, hap1, hap2 ->
            def meta_new = meta + [assembly_stage: "purged"]
            [meta_new, hap1, hap2]
        }

    ch_all_assemblies_after_purging = ch_assemblies_to_purge.no_purge.mix(ch_purged_assemblies)

    //
    // Logic: Filter the input assemblies so that if purging is enabled for an assembly type,
    //        only purged assemblies are polished.
    //
    //        Reshape so that we can concatenate hap1/hap2 together.
    //
    ch_assemblies_to_polish = ch_all_assemblies_after_purging
        | branch { meta, hap1, hap2 ->
            // If we have purged this type of assembly, remove raw stages
            def purging_filter  = (meta.purge == true) ? (meta.assembly_stage != "raw") : meta.assembly_stage == "raw"
            def polish_filter   = (meta.polish == true)
            def polishing_enabled = (params.enable_polishing && params.polishing_longranger_container_path)

            polish: (polishing_enabled && polish_filter && purging_filter)
                return [meta, hap1, hap2]
            no_polish: true
        }

    //
    // Module: Concatenate hap1/hap2 together for polishing
    //
    CONCATENATE_ASSEMBLIES(ch_assemblies_to_polish.polish)
    ch_versions = ch_versions.mix(CONCATENATE_ASSEMBLIES.out.versions)

    //
    // Subworkflow: run polishing
    //
    POLISHING_10X(
        CONCATENATE_ASSEMBLIES.out.file_out,
        ch_i10x_reads
    )
    ch_versions = ch_versions.mix(POLISHING_10X.out.versions)

    //
    // Module: Separate back out primary/alt/hap1/hap2 contigs
    //
    ch_haps = Channel.of("hap1", "hap2")
    ch_assemblies_to_separate = POLISHING_10X.out.assemblies
        | combine(ch_haps)
        | map { meta, asm, hap ->
            [meta + [_hap: hap], asm]
        }

    SEQKIT_GREP_SPLIT_HAPS(
        ch_assemblies_to_separate,
        []
    )
    ch_versions = ch_versions.mix(SEQKIT_GREP_SPLIT_HAPS.out.versions)

    //
    // Logic: Take the split assemblies and re-orgnaise into a [meta, hap1, hap2] format
    //        and add metadata
    //
    ch_polished_assemblies_split = SEQKIT_GREP_SPLIT_HAPS.out.filter
        | map { meta, asm ->
            [meta + [assembly_stage: "polished"], asm]
        }
        | branch { meta, asm ->
            def meta_new = meta - meta.subMap("hap")
            hap1: meta._hap == "hap1"
                return [meta_new, asm]
            hap2: meta._hap == "hap2"
                return [meta_new, asm]
        }

    ch_polished_assemblies = ch_polished_assemblies_split.hap1
        | join(ch_polished_assemblies_split.hap2)

    ch_all_assemblies_after_polishing = ch_assemblies_to_polish.no_polish
        | mix(ch_polished_assemblies)

    //
    // Logic: set up for scaffolding
    //
    ch_assemblies_for_scaffolding_split = ch_assemblies
        | branch { meta, _hap1, _hap2 ->
                def scaffold = false
                if(meta.assembly_stage == "polished") { scaffold = true }
                else if(meta.assembly_stage == "purged" && !meta.polish) { scaffold = true }
                else if(meta.assembly_stage == "raw" && !meta.purge && !meta.polish) { scaffold = true }

                scaffold: params.enable_scaffolding && scaffold
                no_scaffold: true
        }

    //
    // Subworkflow: run hic-mapping and scaffolding
    //
    SCAFFOLDING(
        ch_assemblies_for_scaffolding_split,
        ch_hic_reads,
        params.hic_aligner,
        params.hic_mapping_cram_chunk_size,
        params.cool_bin
    )
    ch_versions = ch_versions.mix(SCAFFOLDING.out.versions)

    ch_all_assemblies_after_scaffolding = ch_assemblies_for_scaffolding_split.scaffold
        | mix(ch_assemblies_for_scaffolding_split.no_scaffold)

    //
    // Subworkflow: calculate assembly QC metrics
    //
    GENOME_STATISTICS(
        ch_all_assemblies_after_scaffolding,
        KMERS.out.fastk,
        KMERS.out.maternal_hapdb,
        KMERS.out.paternal_hapdb,
        busco_lineage,
        busco_lineage_directory
    )
    ch_versions = ch_versions.mix(GENOME_STATISTICS.out.versions)


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
