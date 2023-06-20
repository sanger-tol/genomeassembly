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
    reads_plus_assembly_ch     // [ meta, [reads], [assembly] ], where reads are the pacbio files, and assembly is the primary and alternate asms
    model
    prefix  // [ prefix ] prefix for the output files

    main:
    reads_plus_assembly_ch
        .flatMap { meta, reads, assembly -> reads instanceof List ? reads.collect{ [ meta, reads, assembly ] } : [ [ meta, reads, assembly ] ] }
        .multiMap { meta, reads, assembly -> 
            reads_ch: [ meta, reads ]
            assembly_ch: assembly
        }
        .set { input }

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
    
    if (model) { 
        GET_CALCUTS_PARAMS( model )
        PURGEDUPS_CALCUTS( PURGEDUPS_PBCSTAT.out.stat, GET_CALCUTS_PARAMS.out.cutoffs )
    }
    else {
        PURGEDUPS_CALCUTS( PURGEDUPS_PBCSTAT.out.stat )
    }

    // Split assembly and do self alignment
    reads_plus_assembly_ch
        .map { meta, reads, assembly  -> [ meta, assembly ] }
        .set { minimal_assembly_ch }
    PURGEDUPS_SPLITFA( minimal_assembly_ch )
    minimal_assembly_ch.map{ meta, asm -> Math.ceil(asm.size()/1e9).round() }.set{ idx_num }
    MINIMAP2_ALIGN_ASSEMBLY (
        PURGEDUPS_SPLITFA.out.split_fasta,
        [],    // Trigger read to read alignment
        false, // bam output
        false, // cigar in paf file
        false,  // cigar in bam file
        idx_num
    )
    
    PURGEDUPS_PURGEDUPS(
        PURGEDUPS_PBCSTAT.out.basecov
            .join( PURGEDUPS_CALCUTS.out.cutoff )
            .map { meta, cov, cutoff -> [ meta.findAll { !(it.key in [ 'single_end' ]) }, cov, cutoff ] }
            .join( MINIMAP2_ALIGN_ASSEMBLY.out.paf )
    ) 

    minimal_assembly_ch.join( PURGEDUPS_PURGEDUPS.out.bed )
                       .map { meta, assembly, bed -> prefix ? [[id:meta.id, prefix:prefix], assembly, bed] :
                                                    prefix ? [[id:meta.id], assembly, bed] }
                       .set { ch_getseqs_input }
    PURGEDUPS_GETSEQS( ch_getseqs_input )

    emit:
    pri = PURGEDUPS_GETSEQS.out.purged
    alt = PURGEDUPS_GETSEQS.out.haplotigs
}
