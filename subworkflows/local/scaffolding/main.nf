include { CRAM_MAP_ILLUMINA_HIC as HIC_MAPPING           } from '../../../subworkflows/sanger-tol/cram_map_illumina_hic'
include { BAM_STATS_SAMTOOLS as HIC_MAPPING_STATS        } from '../../../subworkflows/nf-core/bam_stats_samtools/main'
include { FASTA_BAM_SCAFFOLDING_YAHS as SCAFFOLDING_YAHS } from '../../../subworkflows/sanger-tol/fasta_bam_scaffolding_yahs/main'

workflow SCAFFOLDING {
    take:
    ch_assemblies                   // [meta, hap1, hap2]
    val_hic_reads                   // [meta, [reads]]
    val_hic_aligner                 // "bwamem2" or "minimap2"
    val_hic_mapping_cram_chunk_size // int > 1
    val_cool_bin                    // int > 1

    main:
    ch_versions = Channel.empty()

    //
    // Logic: Separate hap1/hap2 into separate channel entries, but tag them
    //
    ch_assemblies_split = ch_assemblies
        | multiMap { meta, hap1, hap2 ->
            hap1: [meta + [_hap: "hap1"], hap1]
            hap2: [meta + [_hap: "hap2"], hap2]
        }

    //
    // Logic: combine Hi-C to assemblies, then spit out new channels
    //        with matching metas
    //
    ch_hic_mapping_inputs = ch_assemblies_split.hap1
        | mix(ch_assemblies_split.hap2)
        | combine(val_hic_reads)
        | multiMap { meta, asm, meta_hic, hic ->
            asm: [ meta, asm ]
            hic: [ meta, hic ]
        }

    //
    // Subworkflow: Map Hi-C data to each assembly
    //
    HIC_MAPPING(
        ch_hic_mapping_inputs.asm,
        ch_hic_mapping_inputs.hic,
        val_hic_aligner,
        val_hic_mapping_cram_chunk_size,
    )
    ch_versions = ch_versions.mix(HIC_MAPPING.out.versions)

    //
    // Subworkflow: Calculate stats for Hi-C mapping
    //
    ch_hic_mapping_stats_input = HIC_MAPPING.out.bam
        | combine(HIC_MAPPING.out.bam.filter { _meta, idx -> idx.getExtension() == "csi" }, by: 0)
        | combine(ch_hic_mapping_inputs.asm, by: 0)
        | multiMap { meta, bam, bai, asm ->
            bam: [ meta, bam, bai ]
            asm: [ meta, asm ]
        }

    HIC_MAPPING_STATS(
        ch_hic_mapping_stats_input.bam,
        ch_hic_mapping_stats_input.asm
    )
    ch_versions = ch_versions.mix(HIC_MAPPING_STATS.out.versions)

    //
    // Subworkflow: scaffold assemblies using yahs and create contact maps
    //
    SCAFFOLDING_YAHS(
        ch_hic_mapping_inputs.asm,
        HIC_MAPPING.out.bam,
        true,
        true,
        true,
        val_cool_bin
    )
    ch_versions   = ch_versions.mix(SCAFFOLDING_YAHS.out.versions)

    //
    // Logic: re-join pairs of assemblies from scaffolding to pass for genome statistics
    //
    ch_assemblies_scaffolded_split = SCAFFOLDING_YAHS.out.scaffolds_fasta
        | branch { meta, assembly ->
            def meta_new = meta - meta.subMap("_hap")
            hap1: meta._hap == "hap1"
                return [ meta_new, assembly ]
            hap2: meta._hap == "hap2"
                return [ meta_new, assembly ]
        }

    ch_assemblies_scaffolded = ch_assemblies_scaffolded_split.hap1
        | join(ch_assemblies_scaffolded_split.hap2)
i
    emit:
    assemblies        = ch_assemblies_scaffolded
    agp               = SCAFFOLDING_YAHS.out.scaffolds_agp
    hic_bin           = SCAFFOLDING_YAHS.out.yahs_bin
    bam               = HIC_MAPPING.out.bam
    yahs_inital       = SCAFFOLDING_YAHS.out.yahs_inital
    yahs_intermediate = SCAFFOLDING_YAHS.out.yahs_intermediate
    yahs_log          = SCAFFOLDING_YAHS.out.yahs_log
    pretext           = SCAFFOLDING_YAHS.out.pretext
    pretext_png       = SCAFFOLDING_YAHS.out.pretext_png
    cool              = SCAFFOLDING_YAHS.out.cool
    hic               = SCAFFOLDING_YAHS.out.hic
    versions   = ch_versions
}
