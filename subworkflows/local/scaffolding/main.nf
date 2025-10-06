include { HIC_MAPPING       } from '../../../subworkflows/sanger-tol/hic_mapping'
include { HIC_MAPPING_STATS } from '../../../subworkflows/local/hic_mapping_stats'
include { SCAFFOLDING_YAHS  } from '../../../subworkflows/local/scaffolding_yahs'

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
        true // mark duplicates
    )
    ch_versions = ch_versions.mix(HIC_MAPPING.out.versions)

    //
    // Subworkflow: Calculate stats for Hi-C mapping
    //
    HIC_MAPPING_STATS(
        HIC_MAPPING.out.bam,
        ch_hic_mapping_inputs.asm
    )
    ch_versions = ch_versions.mix(HIC_MAPPING_STATS.out.versions)

    //
    // Subworkflow: scaffold assemblies using yahs and create contact maps
    //
    SCAFFOLDING_YAHS(
        ch_hic_mapping_inputs.asm,
        HIC_MAPPING.out.bam,
        val_cool_bin
    )
    ch_versions   = ch_versions.mix(SCAFFOLDING_YAHS.out.versions)

    //
    // Logic: re-join pairs of assemblies from scaffolding to pass for genome statistics
    //
    ch_assemblies_scaffolded_split = SCAFFOLDING_YAHS.out.assemblies
        | branch { meta, assembly ->
            def meta_new = meta - meta.subMap("_hap")
            hap1: meta._hap == "hap1"
                return [ meta_new, assembly ]
            hap2: meta._hap == "hap2"
                return [ meta_new, assembly ]
        }

    ch_assemblies_scaffolded = ch_assemblies_scaffolded_split.hap1
        | join(ch_assemblies_scaffolded_split.hap2)

    emit:
    assemblies = ch_assemblies_scaffolded
    bam        = HIC_MAPPING.out.bam
    versions   = ch_versions
}
