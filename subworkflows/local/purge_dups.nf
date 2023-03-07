/*
Adapted from
https://github.com/NBISweden/Earth-Biogenome-Project-pilot/blob/5ec2002638055bb8396857a8ee418bf86188fc59/subworkflows/purge_dups/main.nf
*/

/*
 * Workflow based around the purge_dups tool
 * https://github.com/dfguan/purge_dups
 */

// TODO:: purgedups are no on nf-core
include { MINIMAP2_ALIGN as MINIMAP2_ALIGN_READS    } from "../../modules/nf-core/minimap2/align/main"
include { MINIMAP2_ALIGN as MINIMAP2_ALIGN_ASSEMBLY } from "../../modules/nf-core/minimap2/align/main"
include { PURGEDUPS_CALCUTS                         } from "../../modules/local/purge_dups/calcuts"
include { PURGEDUPS_GETSEQS                         } from "../../modules/local/purge_dups/getseqs"
include { PURGEDUPS_PBCSTAT                         } from "../../modules/local/purge_dups/pbcstat"
include { PURGEDUPS_SPLITFA                         } from "../../modules/local/purge_dups/splitfa"
include { PURGEDUPS_PURGEDUPS                       } from "../../modules/local/purge_dups/purgedups"

workflow PURGE_DUPS {

    take:
    reads_plus_assembly_ch     // [ meta, [reads], [assembly] ], where reads are the pacbio files, and assembly is the primary and alternate asms

    main:
    reads_plus_assembly_ch
        .flatMap { meta, reads, assembly -> reads instanceof List ? reads.collect{ [ meta + [ single_end: true ], it, assembly ] } : [ [ meta + [ single_end: true ], reads, assembly ] ] }
        .multiMap { meta, reads, assembly -> 
            reads_ch: [ meta, reads ]
            assembly_ch: assembly
        }
        .set { input }
    reads_plus_assembly_ch
        .map { meta, reads, assembly -> [ meta, assembly ] }
        .set { assembly_ch }
    // Map pacbio reads
    MINIMAP2_ALIGN_READS(
        input.reads_ch,
        input.assembly_ch,
        false,  // bam output
        false,  // cigar in paf file
        false,  // cigar in bam file
        false  // no split index
    )
    PURGEDUPS_PBCSTAT( MINIMAP2_ALIGN_READS.out.paf.groupTuple() )
    PURGEDUPS_CALCUTS( PURGEDUPS_PBCSTAT.out.stat )
    // TODO:: Cutoffs can likely be estimated from genescope model output.

    // Split assembly and do self alignment
    PURGEDUPS_SPLITFA( assembly_ch )
    assembly_ch.map{ meta, asm -> Math.ceil(asm.length()/1e9).round() }.set{ idx_num }
    MINIMAP2_ALIGN_ASSEMBLY (
        PURGEDUPS_SPLITFA.out.split_fasta,
        [],    // Trigger read to read alignment
        false, // bam output
        false, // cigar in paf file
        false,  // cigar in bam file
        idx_num
    )
    

    // TODO: Check the output from here
    PURGEDUPS_PURGEDUPS(
        PURGEDUPS_PBCSTAT.out.basecov
            .join( PURGEDUPS_CALCUTS.out.cutoff )
            .map { meta, cov, cutoff -> [ meta.findAll { !(it.key in [ 'single_end' ]) }, cov, cutoff ] }
            .join( MINIMAP2_ALIGN_ASSEMBLY.out.paf )
    ) 

    PURGEDUPS_GETSEQS( assembly_ch.join( PURGEDUPS_PURGEDUPS.out.bed ) )

    // TODO: Mix haplotigs back into haplotig set / Verify alternate contigs.

    emit:
    assembly = PURGEDUPS_GETSEQS.out.purged
    coverage = PURGEDUPS_PBCSTAT.out.basecov

}
