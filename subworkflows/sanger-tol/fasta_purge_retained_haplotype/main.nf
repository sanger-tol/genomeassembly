/*
Adapted from
https://github.com/NBISweden/Earth-Biogenome-Project-pilot/blob/5ec2002638055bb8396857a8ee418bf86188fc59/subworkflows/purge_dups/main.nf
*/

include { CAT_CAT as CONCATENATE_HAPLOTYPES         } from '../../../modules/nf-core/cat/cat'
include { MINIMAP2_ALIGN as MINIMAP2_ALIGN_ASSEMBLY } from '../../../modules/nf-core/minimap2/align'
include { PURGEDUPS_CALCUTS                         } from '../../../modules/nf-core/purgedups/calcuts'
include { PURGEDUPS_GETSEQS                         } from '../../../modules/nf-core/purgedups/getseqs'
include { PURGEDUPS_HISTPLOT                        } from '../../../modules/nf-core/purgedups/histplot'
include { PURGEDUPS_PBCSTAT                         } from '../../../modules/nf-core/purgedups/pbcstat'
include { PURGEDUPS_PURGEDUPS                       } from '../../../modules/nf-core/purgedups/purgedups'
include { PURGEDUPS_SPLITFA                         } from '../../../modules/nf-core/purgedups/splitfa'

include { FASTX_MAP_LONG_READS                      } from '../fastx_map_long_reads/main'

workflow FASTA_PURGE_RETAINED_HAPLOTYPE {

    take:
    ch_assemblies             // [meta, hap1, hap2]
    ch_long_reads             // [meta, [reads]]
    val_fastx_reads_per_chunk // integer: number of reads per chunk to map

    main:
    ch_versions = Channel.empty()

    //
    // Logic: split assemblies into primary and alternate
    //
    ch_assemblies_split = ch_assemblies
        .multiMap { meta, hap1, hap2 ->
            primary:  [meta, hap1]
            alternate:[meta, hap2]
        }

    //
    // Subworkflow: Map HiFi reads to primary contigs
    //
    FASTX_MAP_LONG_READS(
        ch_assemblies_split.primary,
        ch_long_reads,
        val_fastx_reads_per_chunk,
        false
    )
    ch_versions = ch_versions.mix(FASTX_MAP_LONG_READS.out.versions)

    //
    // Module: Create read depth histogram
    //
    PURGEDUPS_PBCSTAT(FASTX_MAP_LONG_READS.out.paf)
    ch_versions = ch_versions.mix(PURGEDUPS_PBCSTAT.out.versions)

    //
    // Module: Generate cutoffs based on histogram and kmer coverage
    //
    PURGEDUPS_CALCUTS(PURGEDUPS_PBCSTAT.out.stat)
    ch_versions = ch_versions.mix(PURGEDUPS_CALCUTS.out.versions)

    //
    // Module: Plot purge_dups histogram with cutoffs
    //
    ch_purgedups_histplot_input = PURGEDUPS_PBCSTAT.out.stat
        .combine(PURGEDUPS_CALCUTS.out.cutoff, by: 0)

    PURGEDUPS_HISTPLOT(ch_purgedups_histplot_input)
    ch_versions = ch_versions.mix(PURGEDUPS_HISTPLOT.out.versions)

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
        .combine(PURGEDUPS_CALCUTS.out.cutoff, by: 0)
        .combine(MINIMAP2_ALIGN_ASSEMBLY.out.paf, by: 0)

    PURGEDUPS_PURGEDUPS(ch_purgedups_input)
    ch_versions = ch_versions.mix(PURGEDUPS_PURGEDUPS.out.versions)

    //
    // Module: Generate the primary and alternative contigs
    //
    ch_getseqs_input = ch_assemblies_split.primary
        .combine(PURGEDUPS_PURGEDUPS.out.bed, by: 0)

    PURGEDUPS_GETSEQS(ch_getseqs_input)
    ch_versions = ch_versions.mix(PURGEDUPS_GETSEQS.out.versions)

    //
    // Module: join the alternate assembly to the purged haplotigs
    //         and prepare those that have both for concatenation
    //
    ch_alt_split = ch_assemblies_split.alternate
        .combine(PURGEDUPS_GETSEQS.out.haplotigs, by: 0)
        .branch { meta, haps, alt ->
            def alt_exists  = !!alt
            def haps_exists = haps.size() > 0

            concatenate: alt_exists && haps_exists
                return [ meta, [ alt, haps ] ]
            asis: true
                return [ meta, alt_exists ? alt : haps ]
        }

    //
    // Module: Combine the haplotigs purged from the primary back
    //         into the alternate assembly
    //
    CONCATENATE_HAPLOTYPES(ch_alt_split.concatenate)
    ch_versions = ch_versions.mix(CONCATENATE_HAPLOTYPES.out.versions)

    //
    // Logic: mix the concatenated alts and as-is alts back together
    //
    ch_alts = ch_alt_split.asis.mix(CONCATENATE_HAPLOTYPES.out.file_out)

    emit:
    purged_assemblies          = PURGEDUPS_GETSEQS.out.purged.combine(ch_alts, by: 0)
    purged_haplotigs           = PURGEDUPS_GETSEQS.out.haplotigs
    purgedups_splitfa          = PURGEDUPS_SPLITFA.out.split_fasta
    purgedups_splitfa_self_paf = MINIMAP2_ALIGN_ASSEMBLY.out.paf
    purgedups_pbcstat_hist     = PURGEDUPS_PBCSTAT.out.stat
    purgedups_pbcstat_basecov  = PURGEDUPS_PBCSTAT.out.basecov
    purgedups_calcuts_cutoffs  = PURGEDUPS_CALCUTS.out.cutoff
    purgedups_calcuts_log      = PURGEDUPS_CALCUTS.out.log
    purgedups_histplot         = PURGEDUPS_HISTPLOT.out.png
    purgedups_bed              = PURGEDUPS_PURGEDUPS.out.bed
    purgedups_log              = PURGEDUPS_PURGEDUPS.out.log
    primary_reads_paf          = FASTX_MAP_LONG_READS.out.paf
    versions                   = ch_versions
}
