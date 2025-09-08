/*
Adapted from
https://github.com/NBISweden/Earth-Biogenome-Project-pilot/blob/5ec2002638055bb8396857a8ee418bf86188fc59/subworkflows/purge_dups/main.nf
*/

/*
 * Workflow based around the purge_dups tool
 * https://github.com/dfguan/purge_dups
 */

include { CAT_CAT as CAT_PURGED_HAPS_TO_ALT         } from '../../../modules/nf-core/cat/cat'
include { MINIMAP2_ALIGN as MINIMAP2_ALIGN_READS    } from '../../../modules/nf-core/minimap2/align'
include { MINIMAP2_ALIGN as MINIMAP2_ALIGN_ASSEMBLY } from '../../../modules/nf-core/minimap2/align'
include { PURGEDUPS_CALCUTS                         } from '../../../modules/nf-core/purgedups/calcuts'
include { PURGEDUPS_GETSEQS                         } from '../../../modules/nf-core/purgedups/getseqs'
include { PURGEDUPS_PBCSTAT                         } from '../../../modules/nf-core/purgedups/pbcstat'
include { PURGEDUPS_PURGEDUPS                       } from '../../../modules/nf-core/purgedups/purgedups'
include { PURGEDUPS_SPLITFA                         } from '../../../modules/nf-core/purgedups/splitfa'

workflow PURGING {

    take:
    ch_assemblies     // [meta, hap1, hap2]
    val_long_reads    // [meta, [reads]] - should be a value channel

    main:
    ch_versions = Channel.empty()

    //
    // Logic: split assemblies into primary and alternate
    //

    ch_assemblies_split = ch_assemblies
        | multiMap { meta, hap1, hap2 ->
            primary:  [meta, hap1]
            alternate:[meta, hap2]
        }

    //
    // Module: Map HiFi reads to primary contigs
    //
    ch_input_for_minimap2_align_reads = ch_assemblies_split.primary
        | combine(val_long_reads)
        | multiMap { meta, asm, meta_reads, reads ->
            reads: [meta, reads]
            asm:   [meta, asm]
        }

    MINIMAP2_ALIGN_READS(
        ch_input_for_minimap2_align_reads.reads,
        ch_input_for_minimap2_align_reads.asm,
        false,      // bam output
        false,      // bam index extension
        false,      // cigar in paf file
        false       // cigar in bam file
    )
    ch_versions = ch_versions.mix(MINIMAP2_ALIGN_READS.out.versions)

    //
    // Module: Create read depth histogram
    //
    PURGEDUPS_PBCSTAT(MINIMAP2_ALIGN_READS.out.paf)
    ch_versions = ch_versions.mix(PURGEDUPS_PBCSTAT.out.versions)

    //
    // Module: Generate cutoffs based on histogram and kmer coverage
    //
    PURGEDUPS_CALCUTS(PURGEDUPS_PBCSTAT.out.stat)
    ch_versions = ch_versions.mix(PURGEDUPS_CALCUTS.out.versions)

    //
    // Module: Split assembly
    //
    PURGEDUPS_SPLITFA(ch_assemblies_split.primary)
    ch_versions = ch_versions.mix(PURGEDUPS_SPLITFA.out.versions)

    //
    // MODULE: PEFORM SELF ALIGNMENT
    //
    MINIMAP2_ALIGN_ASSEMBLY (
        PURGEDUPS_SPLITFA.out.split_fasta,
        [[:], []], // Trigger read to read alignment
        false,     // bam output
        false,     // bam index extension
        false,     // cigar in paf file
        false      // cigar in bam file
    )
    ch_versions = ch_versions.mix(MINIMAP2_ALIGN_ASSEMBLY.out.versions)

    //
    // Module: Purge haplotigs from primary assembly
    //
    ch_purgedups_input = PURGEDUPS_PBCSTAT.out.basecov
        | join(PURGEDUPS_CALCUTS.out.cutoff)
        | join(MINIMAP2_ALIGN_ASSEMBLY.out.paf)

    PURGEDUPS_PURGEDUPS(ch_purgedups_input)
    ch_versions = ch_versions.mix(PURGEDUPS_PURGEDUPS.out.versions)

    //
    // Module: Generate the primary and alternative contigs
    //
    ch_getseqs_input = ch_assemblies_split.primary
        | join(PURGEDUPS_PURGEDUPS.out.bed)

    PURGEDUPS_GETSEQS(ch_getseqs_input)
    ch_versions = ch_versions.mix(PURGEDUPS_GETSEQS.out.versions)

    //
    // Module: Combine the haplotigs purged from the primary back
    //         into the alternate assembly
    //
    ch_alt_to_cat = PURGEDUPS_GETSEQS.out.haplotigs
        | join(ch_assemblies_split.alternate)
        | map { meta, haps, alt -> [meta, [alt, haps]] }

    CAT_PURGED_HAPS_TO_ALT(ch_alt_to_cat)
    ch_versions = ch_versions.mix(CAT_PURGED_HAPS_TO_ALT.out.versions)

    //
    // Logic: Combine purged primary and alternate assemblies back into
    //        a single output channel
    //
    ch_assemblies  = PURGEDUPS_GETSEQS.out.purged
        | join(CAT_PURGED_HAPS_TO_ALT.out.file_out)

    emit:
    assemblies = ch_assemblies
    versions = ch_versions
}
