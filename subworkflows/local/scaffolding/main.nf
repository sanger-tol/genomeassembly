include { CRAM_MAP_ILLUMINA_HIC                 } from '../../../subworkflows/sanger-tol/cram_map_illumina_hic'
include { BAM_STATS_SAMTOOLS                    } from '../../../subworkflows/nf-core/bam_stats_samtools'
include { FASTA_BAM_SCAFFOLDING_YAHS            } from '../../../subworkflows/sanger-tol/fasta_bam_scaffolding_yahs'
include { PAIRS_CREATE_CONTACT_MAPS             } from '../../../subworkflows/sanger-tol/pairs_create_contact_maps/main'

include { HTSLIB_BGZIPTABIX as BGZIP_SCAFFOLDED } from '../../../modules/nf-core/htslib/bgziptabix'

workflow SCAFFOLDING {
    take:
    ch_scaffolding_specs // spec
    ch_assemblies // [meta, hap1, hap2]
    val_hic_aligner // "bwamem2" or "minimap2"
    val_hic_mapping_cram_chunk_size // int > 1
    val_cool_bin // int > 1
    val_build_pretext_map
    val_build_juicer_map
    val_build_cooler_map

    main:
    //
    // Logic: join all the assemblies with the scaffolding specifications and
    // data, filter for those assemblies which are to be scaffolded.
    //
    // Here we keep the original spec for Hi-C mapping. This ensures that any changes
    // to the scaffolding parameters don't re-run Hi-C mapping on resume.
    //
    ch_hic_mapping_inputs = ch_assemblies
        .combine(ch_scaffolding_specs)
        .filter { input_spec, _asm1, _asm2, spec -> input_spec.id == spec.prevID }
        .multiMap { input_spec, asm1, asm2, spec ->
            def spec_hap1 = input_spec + [_hap: "hap1"]
            def spec_hap2 = input_spec + [_hap: "hap2"]
            hap1: [spec_hap1, asm1]
            hap2: [spec_hap2, asm2]
            hic_reads: [[spec_hap1, spec_hap2], spec.data.hic.reads]
        }

    //
    // Subworkflow: Map Hi-C data to each assembly
    //
    CRAM_MAP_ILLUMINA_HIC(
        ch_hic_mapping_inputs.hap1.mix(ch_hic_mapping_inputs.hap2),
        ch_hic_mapping_inputs.hic_reads.transpose(by: 0),
        val_hic_aligner,
        val_hic_mapping_cram_chunk_size,
    )

    //
    // Logic: Once Hi-C mapping has run, combine everything with the scaffolding specifications
    // so we can attach the scaffolding parameters, and make joining below cleaner by doing this
    // once.
    //
    ch_post_hic_inputs = CRAM_MAP_ILLUMINA_HIC.out.bam
        .combine(CRAM_MAP_ILLUMINA_HIC.out.bam_index.filter { _meta, idx -> idx.getExtension() == "csi" }, by: 0)
        .combine(ch_hic_mapping_inputs.hap1.mix(ch_hic_mapping_inputs.hap2), by: 0)
        .combine(ch_scaffolding_specs)
        .filter { input_spec, _bam, _csi, hap, scaf_spec -> input_spec.id == scaf_spec.prevID }
        .multiMap { input_spec, bam, csi, hap, spec ->
            def out_spec = spec + input_spec.subMap("_hap")
            asm: [out_spec, hap]
            asm_fai: [out_spec, hap, []]
            bam: [out_spec, bam]
            bam_csi: [out_spec, bam, csi]
        }

    //
    // Subworkflow: Calculate stats for Hi-C mapping
    //
    BAM_STATS_SAMTOOLS(
        ch_post_hic_inputs.bam_csi,
        ch_post_hic_inputs.asm_fai
    )

    //
    // Subworkflow: scaffold assemblies using yahs and create contact maps
    //
    // Here we take the Hi-C mapping outputs and replace the spec to
    // attach the scaffolding parameters.
    //
    FASTA_BAM_SCAFFOLDING_YAHS(
        ch_post_hic_inputs.asm,
        ch_post_hic_inputs.bam
    )

    //
    // Subworkflow: create contact maps
    //
    ch_contact_map_inputs = FASTA_BAM_SCAFFOLDING_YAHS.out.scaffolds_pairs
        .join(FASTA_BAM_SCAFFOLDING_YAHS.out.scaffolds_chromsizes)
        .multiMap { meta, pairs, sizes ->
            pairs: [meta, pairs]
            sizes: [meta, sizes]
        }

    PAIRS_CREATE_CONTACT_MAPS(
        ch_contact_map_inputs.pairs,
        ch_contact_map_inputs.sizes,
        channel.empty(),
        val_build_pretext_map,
        val_build_pretext_map ?: false,
        val_build_cooler_map,
        val_build_juicer_map,
        val_cool_bin,
    )

    //
    // Module: bgzip all scaffolded assembly fasta
    //
    BGZIP_SCAFFOLDED(
        FASTA_BAM_SCAFFOLDING_YAHS.out.scaffolds_fasta.map { meta, fasta -> [meta, fasta, [], []] },
        "compress",
        false,
        "fa",
    )

    //
    // Logic: re-join pairs of assemblies from scaffolding to pass for genome statistics
    //
    ch_assemblies_scaffolded = FASTA_BAM_SCAFFOLDING_YAHS.out.scaffolds_fasta
        .filter { meta, _scaffolds -> meta._hap == "hap1" }
        .mix(FASTA_BAM_SCAFFOLDING_YAHS.out.scaffolds_fasta.filter { meta, _scaffolds -> meta._hap == "hap2" })
        .map { meta, asm -> [meta - meta.subMap("_hap"), asm] }
        .groupTuple(size: 2)
        .map { meta, asms -> [meta, asms[0], asms[1]] }

    //
    // Logic: combine all scaffolding outputs into a single map for ease of publishing
    //
    ch_scaffolding_output = BGZIP_SCAFFOLDED.out.output
        .join(ch_post_hic_inputs.bam_csi, by: 0)
        .join(BAM_STATS_SAMTOOLS.out.stats, by: 0)
        .join(BAM_STATS_SAMTOOLS.out.flagstat, by: 0)
        .join(BAM_STATS_SAMTOOLS.out.idxstats, by: 0)
        .join(FASTA_BAM_SCAFFOLDING_YAHS.out.scaffolds_agp, by: 0)
        .join(FASTA_BAM_SCAFFOLDING_YAHS.out.yahs_bin, by: 0)
        .join(FASTA_BAM_SCAFFOLDING_YAHS.out.yahs_inital, by: 0, remainder: true)
        .join(FASTA_BAM_SCAFFOLDING_YAHS.out.yahs_intermediate, by: 0, remainder: true)
        .join(FASTA_BAM_SCAFFOLDING_YAHS.out.yahs_log, by: 0)
        .join(PAIRS_CREATE_CONTACT_MAPS.out.pretext, by: 0, remainder: true)
        .join(PAIRS_CREATE_CONTACT_MAPS.out.pretext_png, by: 0, remainder: true)
        .join(PAIRS_CREATE_CONTACT_MAPS.out.cool, by: 0, remainder: true)
        .join(PAIRS_CREATE_CONTACT_MAPS.out.hic, by: 0, remainder: true)
        .map { spec, fasta, bam, bai, stats, flagstats, idxstats, agp, bin, initial, intermed, log, pretext, png, cool, hic ->
            return spec.subMap(["id", "stage", "data", "params", "tools"]) + [hap: spec._hap, output: [scaffolding: [fasta: fasta, bam: bam, bai: bai, stats: stats, flagstats: flagstats, idxstats: idxstats, yahs_agp: agp, yahs_bin: bin, yahs_initial: initial, yahs_intermeriate: intermed, yahs_log: log, pretext: pretext, pretext_png: png, cool: cool, hic: hic]]]
        }

    emit:
    scaffolded_assemblies = ch_assemblies_scaffolded
    scaffolding_output    = ch_scaffolding_output
}
