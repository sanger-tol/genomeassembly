/*
Adapted from
https://github.com/NBISweden/Earth-Biogenome-Project-pilot/blob/5ec2002638055bb8396857a8ee418bf86188fc59/subworkflows/purge_dups/main.nf
*/

/*
 * Workflow based around the purge_dups tool
 * https://github.com/dfguan/purge_dups
 */

include { MINIMAP2_ALIGN as MINIMAP2_ALIGN_READS    } from "../../modules/nf-core/minimap2/align/main"
include { MINIMAP2_ALIGN as MINIMAP2_ALIGN_ASSEMBLY } from "../../modules/nf-core/minimap2/align/main"

include { PURGEDUPS_CALCUTS                         } from '../../modules/nf-core/purgedups/calcuts/main'
include { PURGEDUPS_GETSEQS                         } from '../../modules/nf-core/purgedups/getseqs/main'
include { PURGEDUPS_PBCSTAT                         } from '../../modules/nf-core/purgedups/pbcstat/main'
include { PURGEDUPS_PURGEDUPS                       } from '../../modules/nf-core/purgedups/purgedups/main'
include { PURGEDUPS_SPLITFA                         } from '../../modules/nf-core/purgedups/splitfa/main'

include { GET_CALCUTS_PARAMS                        } from "../../modules/local/get_calcuts_params"

workflow PURGE_DUPS {

    take:
    reads_plus_assembly_ch     // [ meta, [reads], [assembly], [model] ], where reads are the pacbio files, and assembly is the primary and alternate asms
    prefix  // [ prefix ] prefix for the output files

    main:
    ch_versions = Channel.empty()

    //
    // LOGIC: TRANSFROM INPUT DATA STRUCTURE
    //
    reads_plus_assembly_ch
        .flatMap { meta, reads, assembly, model -> reads instanceof List ? reads.collect{ [ meta, reads, assembly, model ] } : [ [ meta, reads, assembly, model ] ] }
        .multiMap { meta, reads, assembly, model -> 
            reads_ch: [ meta, reads ]
            assembly_ch: assembly
            model_ch: [ meta, model ] 
        }
        .set { input }

    //
    // MODULE: MAP HIFI READS TO CONTIGS
    //
    MINIMAP2_ALIGN_READS(
        input.reads_ch,
        input.assembly_ch,
        false,  // bam output
        false,  // cigar in paf file
        false,  // cigar in bam file
        false  // no split index
    )
    ch_versions = ch_versions.mix(MINIMAP2_ALIGN_READS.out.versions)

    //
    // MODULE: CREATE READ DEPTH HISTOGRAM 
    //
    PURGEDUPS_PBCSTAT( MINIMAP2_ALIGN_READS.out.paf.groupTuple() )
    ch_versions = ch_versions.mix(PURGEDUPS_PBCSTAT.out.versions)
    
    //
    // MODULE: PARSE KMER COVERAGE  
    //
    GET_CALCUTS_PARAMS( input.model_ch )
    ch_versions = ch_versions.mix(GET_CALCUTS_PARAMS.out.versions)

    //
    // MODULE: GENERATES CUTOFFS BASED ON HISTOGRAM AND KMER COVERAGE
    //
    PURGEDUPS_CALCUTS( PURGEDUPS_PBCSTAT.out.stat, GET_CALCUTS_PARAMS.out.cutoffs )
    ch_versions = ch_versions.mix(PURGEDUPS_CALCUTS.out.versions)

    //
    // LOGIC: TRANSFORM ASSEMBLY INPUT
    //
    reads_plus_assembly_ch
        .map { meta, reads, assembly, model -> [ meta, assembly ] }
        .set { minimal_assembly_ch }

    //
    // MODULE: SPLIT ASSEMLBY
    //
    PURGEDUPS_SPLITFA( minimal_assembly_ch )
    ch_versions = ch_versions.mix(PURGEDUPS_SPLITFA.out.versions)

    //
    // LOGIC: ESTIMATE THE SIZE OF THE INDEX
    //
    minimal_assembly_ch.map{ meta, asm -> Math.ceil(asm.size()/1e9).round() }.set{ idx_num }

    //
    // MODULE: PEFORM SELF ALIGNMENT
    //
    MINIMAP2_ALIGN_ASSEMBLY (
        PURGEDUPS_SPLITFA.out.split_fasta,
        [],    // Trigger read to read alignment
        false, // bam output
        false, // cigar in paf file
        false,  // cigar in bam file
        idx_num
    )
    ch_versions = ch_versions.mix(MINIMAP2_ALIGN_ASSEMBLY.out.versions)
    
    //
    // MODULE: PURGE HAPLOTIGS
    //
    PURGEDUPS_PURGEDUPS(
        PURGEDUPS_PBCSTAT.out.basecov
            .join( PURGEDUPS_CALCUTS.out.cutoff )
            .map { meta, cov, cutoff -> [ meta.findAll { !(it.key in [ 'single_end' ]) }, cov, cutoff ] }
            .join( MINIMAP2_ALIGN_ASSEMBLY.out.paf )
    ) 
    ch_versions = ch_versions.mix(PURGEDUPS_PURGEDUPS.out.versions)

    //
    // LOGIC: PREPARE INPUT FOR GETSEQS 
    //
    minimal_assembly_ch.join( PURGEDUPS_PURGEDUPS.out.bed )
                       .map { meta, assembly, bed -> [[id:meta.id, prefix:prefix], assembly, bed] }
                       .set { ch_getseqs_input }

    //
    // MODULE: GENERATE PRIMARY AND ALT CONTIGS
    //
    PURGEDUPS_GETSEQS( ch_getseqs_input )
    ch_versions = ch_versions.mix(PURGEDUPS_GETSEQS.out.versions)

    emit:
    pri = PURGEDUPS_GETSEQS.out.purged
    alt = PURGEDUPS_GETSEQS.out.haplotigs

    versions = ch_versions
}
