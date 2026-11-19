/*
Adapted from
https://github.com/NBISweden/Earth-Biogenome-Project-pilot/blob/5ec2002638055bb8396857a8ee418bf86188fc59/subworkflows/purge_dups/main.nf
*/

include { FIND_CONCATENATE as CONCATENATE_HAPLOTYPES } from '../../../modules/nf-core/find/concatenate'
include { MINIMAP2_ALIGN as MINIMAP2_ALIGN_ASSEMBLY  } from '../../../modules/nf-core/minimap2/align'
include { PURGEDUPS_CALCUTS                          } from '../../../modules/nf-core/purgedups/calcuts'
include { PURGEDUPS_GETSEQS                          } from '../../../modules/nf-core/purgedups/getseqs'
include { PURGEDUPS_HISTPLOT                         } from '../../../modules/nf-core/purgedups/histplot'
include { PURGEDUPS_PBCSTAT                          } from '../../../modules/nf-core/purgedups/pbcstat'
include { PURGEDUPS_PURGEDUPS                        } from '../../../modules/nf-core/purgedups/purgedups'
include { PURGEDUPS_SPLITFA                          } from '../../../modules/nf-core/purgedups/splitfa'

include { FASTX_MAP_LONG_READS                       } from '../fastx_map_long_reads/main'

workflow FASTA_PURGE_RETAINED_HAPLOTYPE {
    take:
    ch_assemblies // [meta, [hap1, hap2], ..., hap_n]
    ch_long_reads // [meta, [reads]]
    val_fastx_reads_per_chunk // integer: number of reads per chunk to map

    main:
    //
    // Logic: validate that input assemblies are unzipped
    //
    ch_assemblies.subscribe { _meta, haps ->
        if (haps.collect { hap -> hap.getExtension() == "gz" }.any()) {
            error("Error: Input assemblies to FASTA_PURGE_RETAINED_HAPLOTYPE must be unzipped!")
        }
    }

    //
    // Logic: split assemblies into primary and alternate
    //
    ch_assemblies_split = ch_assemblies.multiMap { meta, haps ->
        primary: [meta, haps[0]]
        alternate: [meta, haps.drop(1)]
    }

    //
    // Subworkflow: Map HiFi reads to primary contigs
    //
    FASTX_MAP_LONG_READS(
        ch_assemblies_split.primary,
        ch_long_reads,
        val_fastx_reads_per_chunk,
        false,
        channel.empty(),
    )

    //
    // Module: Create read depth histogram
    //
    PURGEDUPS_PBCSTAT(FASTX_MAP_LONG_READS.out.paf)

    //
    // Module: Generate cutoffs based on histogram and kmer coverage
    //
    PURGEDUPS_CALCUTS(PURGEDUPS_PBCSTAT.out.stat)

    //
    // Module: Plot purge_dups histogram with cutoffs
    //
    ch_purgedups_histplot_input = PURGEDUPS_PBCSTAT.out.stat.combine(PURGEDUPS_CALCUTS.out.cutoff, by: 0)

    PURGEDUPS_HISTPLOT(ch_purgedups_histplot_input)

    //
    // Module: Split assembly
    //
    PURGEDUPS_SPLITFA(ch_assemblies_split.primary)

    //
    // MODULE: PEFORM SELF ALIGNMENT
    //
    MINIMAP2_ALIGN_ASSEMBLY(
        PURGEDUPS_SPLITFA.out.split_fasta,
        [[:], []], // Trigger read to read alignment
        false,     // bam output
        false,     // bam index extension
        false,     // cigar in paf file
        false      // cigar in bam file
    )

    //
    // Module: Purge haplotigs from primary assembly
    //
    ch_purgedups_input = PURGEDUPS_PBCSTAT.out.basecov
        .combine(PURGEDUPS_CALCUTS.out.cutoff, by: 0)
        .combine(MINIMAP2_ALIGN_ASSEMBLY.out.paf, by: 0)

    PURGEDUPS_PURGEDUPS(ch_purgedups_input)

    //
    // Module: Generate the primary and alternative contigs
    //
    ch_getseqs_input = ch_assemblies_split.primary.combine(PURGEDUPS_PURGEDUPS.out.bed, by: 0)

    PURGEDUPS_GETSEQS(ch_getseqs_input)

    //
    // Module: join the alternate assembly to the purged haplotigs
    //         and prepare those that have both for concatenation
    //
    ch_alt_split = ch_assemblies_split.alternate
        .join(PURGEDUPS_GETSEQS.out.haplotigs, by: 0, remainder: true)
        .map { meta, input_alts, purged_haplotigs -> [meta, (input_alts << purged_haplotigs).findAll { it -> it.size() > 0 }] }
        .branch { _meta, alts ->
            concatenate: alts.size() > 1
            asis: true
        }

    //
    // Module: Combine the haplotigs purged from the primary back
    //         into the alternate assembly
    //
    CONCATENATE_HAPLOTYPES(ch_alt_split.concatenate)

    //
    // Logic: mix the concatenated alts and as-is alts back together
    //
    ch_alts = ch_alt_split.asis.mix(CONCATENATE_HAPLOTYPES.out.file_out)

    //
    // Logic: combine the purged primary and alternate assemblies into a single output channel
    //
    ch_purged_assembly_set = PURGEDUPS_GETSEQS.out.purged
        .join(ch_alts, by: 0, remainder: true)
        .map { meta, purged, alts ->
            [meta, [purged, alts].flatten().findAll { it -> it.size() > 0 }]
        }

    emit:
    purged_assemblies          = ch_purged_assembly_set
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
}
