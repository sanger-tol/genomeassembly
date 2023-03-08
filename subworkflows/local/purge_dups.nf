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
include { GET_CALCUTS_PARAMS                        } from "../../modules/local/get_calcuts_params"

include { CAT_CAT                                   } from "../../modules/nf-core/cat/cat/main"
include { FASTK_FASTK                               } from "../../modules/nf-core/fastk/fastk/main"
include { FASTK_HISTEX                              } from '../../modules/nf-core/fastk/histex/main'
include { GENESCOPEFK                               } from "../../modules/nf-core/genescopefk/main"

workflow PURGE_DUPS {

    take:
    reads_plus_assembly_ch     // [ meta, [reads], [assembly] ], where reads are the pacbio files, and assembly is the primary and alternate asms

    main:
    reads_plus_assembly_ch
        .flatMap { meta, reads, assembly -> reads instanceof List ? reads.collect{ [ meta + [ single_end: true], it, assembly ] } : [ [ meta + [ single_end: true], reads, assembly ] ] }
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
    CAT_CAT( input.reads_ch )
    CAT_CAT.out.file_out.map{ meta, path -> path.renameTo(path.getParent() + '/' + path.getBaseName() + '.fa.gz'); 
                                            [meta, path.getParent() + '/' + path.getBaseName() + '.fa.gz'] }
                        .set{ reads_merged_ch }
    FASTK_FASTK( reads_merged_ch )
    FASTK_HISTEX( FASTK_FASTK.out.hist )
    GENESCOPEFK ( FASTK_HISTEX.out.hist )
    GET_CALCUTS_PARAMS( GENESCOPEFK.out.model )
    GET_CALCUTS_PARAMS.out.cutoffs.view()

    PURGEDUPS_CALCUTS( PURGEDUPS_PBCSTAT.out.stat, GET_CALCUTS_PARAMS.out.cutoffs )

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
