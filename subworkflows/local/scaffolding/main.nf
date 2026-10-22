include { CRAM_MAP_ILLUMINA_HIC                 } from '../../../subworkflows/sanger-tol/cram_map_illumina_hic'
include { BAM_STATS_SAMTOOLS                    } from '../../../subworkflows/nf-core/bam_stats_samtools'
include { FASTA_BAM_SCAFFOLDING_YAHS            } from '../../../subworkflows/sanger-tol/fasta_bam_scaffolding_yahs'
include { FASTA_BAM_SCAFFOLDING_RAGTAG_YAHS     } from '../../../subworkflows/sanger-tol/fasta_bam_scaffolding_ragtag_yahs/main'
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
        .filter { input_spec, _bam, _csi, _hap, scaf_spec -> input_spec.id == scaf_spec.prevID }
        .map { input_spec, bam, csi, hap, scaf_spec ->
            [scaf_spec + input_spec.subMap("_hap"), hap, bam, csi]
        }

    //
    // Subworkflow: Calculate stats for Hi-C mapping
    //
    ch_bam_stats_samtools_input = ch_post_hic_inputs.multiMap { spec, hap, bam, csi ->
        bam_csi: [spec, bam, csi]
        asm_fai: [spec, hap, []]
    }
    BAM_STATS_SAMTOOLS(
        ch_bam_stats_samtools_input.bam_csi,
        ch_bam_stats_samtools_input.asm_fai,
    )

    //
    // Subworkflow: scaffold assemblies using yahs and create contact maps.
    //
    // We pass in all assemblies to be scaffolded in "yahs" mode, and the hap1
    // assembly if the mode is "ragtag_yahs" AND there is not a pre-provided
    // RagTag reference.
    //
    ch_fasta_bam_scaffolding_yahs_inputs = ch_post_hic_inputs
        .filter { spec, _hap, _bam, _csi ->
            spec.params.scaffolding_mode == "yahs" || [
                spec.params.scaffolding_mode == "ragtag_yahs",
                spec._hap == "hap1",
                !spec.params.ratag_reference,
            ].every()
        }
        .multiMap { spec, hap, bam, _csi ->
            asm: [spec, hap]
            bam: [spec, bam]
        }

    FASTA_BAM_SCAFFOLDING_YAHS(
        ch_fasta_bam_scaffolding_yahs_inputs.asm,
        ch_fasta_bam_scaffolding_yahs_inputs.bam,
    )

    //
    // Logic: identify the RagTag reference. One of the two of these MUST always exist:
    //     - a pre-provided reference in spec.params.ratag_reference
    //     - the YaHS-scaffolded hap1 assembly from the previous step
    //
    // So mixing the two will always produce a channel of 1 for each assembly specification.
    //
    ch_ragtag_reference = ch_post_hic_inputs
        .filter { spec, _hap, _bam, _csi ->
            spec.params.scaffolding_mode == "ragtag_yahs" && spec._hap == "hap1" && spec.params.ratag_reference
        }
        .map { spec, _hap, _bam, _csi -> [spec - spec.subMap("_hap"), spec.params.ratag_reference] }
        .mix(
            FASTA_BAM_SCAFFOLDING_YAHS.out.scaffolds_fasta.filter { spec, _fasta ->
                spec.params.scaffolding_mode == "ragtag_yahs" && spec._hap == "hap1" && !spec.params.ratag_reference
            }.map { spec, fasta -> [spec - spec.subMap("_hap"), fasta] }
        )

    //
    // Subworkflow: scaffold assemblies using ragtag and yahs.
    //
    // Generate inputs by combining the assemblies with the RagTag reference.
    // We have to be careful to filter out the hap1 assembly if it was already
    // scaffolded with YaHS.
    //
    ch_fasta_bam_scaffolding_ragtag_yahs_inputs = ch_post_hic_inputs
        .filter { spec, _hap, _bam, _csi ->
            [
                spec.params.scaffolding_mode == "ragtag_yahs",
                spec.params.ratag_reference ? true : spec._hap != "hap1",
            ].every()
        }
        .combine(ch_ragtag_reference)
        .filter { spec, _hap, _bam, _csi, ref_spec, _ref_asm -> spec.id == ref_spec.id }
        .multiMap { spec, hap, bam, _csi, _ref_spec, reference ->
            asm: [spec, hap]
            bam: [spec, bam]
            ref: [spec, reference]
        }

    FASTA_BAM_SCAFFOLDING_RAGTAG_YAHS(
        ch_fasta_bam_scaffolding_ragtag_yahs_inputs.asm,
        ch_fasta_bam_scaffolding_ragtag_yahs_inputs.bam,
        ch_fasta_bam_scaffolding_ragtag_yahs_inputs.ref,
    )

    //
    // Logic: Combine all YaHS outputs into mixed channels so we can forget whether we ran
    // YaHS or RagTag+YaHS scaffolding.
    //
    ch_yahs_scaffolds = FASTA_BAM_SCAFFOLDING_YAHS.out.scaffolds_fasta.mix(FASTA_BAM_SCAFFOLDING_RAGTAG_YAHS.out.scaffolds_fasta)
    ch_yahs_scaffolds_agp = FASTA_BAM_SCAFFOLDING_YAHS.out.scaffolds_agp.mix(FASTA_BAM_SCAFFOLDING_RAGTAG_YAHS.out.scaffolds_agp)
    ch_yahs_scaffolds_pairs = FASTA_BAM_SCAFFOLDING_YAHS.out.scaffolds_pairs.mix(FASTA_BAM_SCAFFOLDING_RAGTAG_YAHS.out.scaffolds_pairs)
    ch_yahs_scaffolds_chromsizes = FASTA_BAM_SCAFFOLDING_YAHS.out.scaffolds_chromsizes.mix(FASTA_BAM_SCAFFOLDING_RAGTAG_YAHS.out.scaffolds_chromsizes)
    ch_yahs_bin = FASTA_BAM_SCAFFOLDING_YAHS.out.yahs_bin.mix(FASTA_BAM_SCAFFOLDING_RAGTAG_YAHS.out.yahs_bin)
    ch_yahs_inital = FASTA_BAM_SCAFFOLDING_YAHS.out.yahs_inital.mix(FASTA_BAM_SCAFFOLDING_RAGTAG_YAHS.out.yahs_inital)
    ch_yahs_intermediate = FASTA_BAM_SCAFFOLDING_YAHS.out.yahs_intermediate.mix(FASTA_BAM_SCAFFOLDING_RAGTAG_YAHS.out.yahs_intermediate)
    ch_yahs_log = FASTA_BAM_SCAFFOLDING_YAHS.out.yahs_log.mix(FASTA_BAM_SCAFFOLDING_RAGTAG_YAHS.out.yahs_log)

    //
    // Subworkflow: create contact maps
    //
    ch_contact_map_inputs = ch_yahs_scaffolds_pairs
        .join(ch_yahs_scaffolds_chromsizes)
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
        ch_yahs_scaffolds.map { meta, fasta -> [meta, fasta, [], []] },
        "compress",
        false,
        "fa",
    )

    //
    // Logic: re-join pairs of assemblies from scaffolding to pass for genome statistics
    //
    ch_assemblies_scaffolded = ch_yahs_scaffolds
        .filter { meta, _scaffolds -> meta._hap == "hap1" }
        .mix(FASTA_BAM_SCAFFOLDING_YAHS.out.scaffolds_fasta.filter { meta, _scaffolds -> meta._hap == "hap2" })
        .map { meta, asm -> [meta - meta.subMap("_hap"), asm] }
        .groupTuple(size: 2)
        .map { meta, asms -> [meta, asms[0], asms[1]] }

    //
    // Logic: combine all scaffolding outputs into a single map for ease of publishing
    //
    ch_scaffolding_output = BGZIP_SCAFFOLDED.out.output
        .join(ch_bam_stats_samtools_input.bam_csi, by: 0)
        .join(BAM_STATS_SAMTOOLS.out.stats, by: 0)
        .join(BAM_STATS_SAMTOOLS.out.flagstat, by: 0)
        .join(BAM_STATS_SAMTOOLS.out.idxstats, by: 0)
        .join(ch_yahs_scaffolds_agp, by: 0)
        .join(ch_yahs_bin, by: 0)
        .join(ch_yahs_inital, by: 0, remainder: true)
        .join(ch_yahs_intermediate, by: 0, remainder: true)
        .join(ch_yahs_log, by: 0)
        .join(FASTA_BAM_SCAFFOLDING_RAGTAG_YAHS.out.ragtag_scaffolds, by: 0, remainder: true)
        .join(FASTA_BAM_SCAFFOLDING_RAGTAG_YAHS.out.ragtag_agp, by: 0, remainder: true)
        .join(FASTA_BAM_SCAFFOLDING_RAGTAG_YAHS.out.ragtag_stats, by: 0, remainder: true)
        .join(FASTA_BAM_SCAFFOLDING_RAGTAG_YAHS.out.ragtag_paf, by: 0, remainder: true)
        .join(FASTA_BAM_SCAFFOLDING_RAGTAG_YAHS.out.ragtag_delta, by: 0, remainder: true)
        .join(FASTA_BAM_SCAFFOLDING_RAGTAG_YAHS.out.confidence_txt, by: 0, remainder: true)
        .join(PAIRS_CREATE_CONTACT_MAPS.out.pretext, by: 0, remainder: true)
        .join(PAIRS_CREATE_CONTACT_MAPS.out.pretext_png, by: 0, remainder: true)
        .join(PAIRS_CREATE_CONTACT_MAPS.out.cool, by: 0, remainder: true)
        .join(PAIRS_CREATE_CONTACT_MAPS.out.hic, by: 0, remainder: true)
        .map { spec, fasta, bam, bai, stats, flagstats, idxstats, yahs_agp, yahs_bin, yahs_initial, yahs_intermed, yahs_log, rt_scaf, rt_agp, rt_stats, rt_paf, rt_delta, rt_conf, pretext, png, cool, hic ->
            return spec.subMap(["id", "stage", "data", "params", "tools"]) + [hap: spec._hap, output: [scaffolding: [
                fasta: fasta,
                bam: bam,
                bai: bai,
                stats: stats,
                flagstats: flagstats,
                idxstats: idxstats,
                yahs_agp: yahs_agp,
                yahs_bin: yahs_bin,
                yahs_initial: yahs_initial,
                yahs_intermeriate: yahs_intermed,
                yahs_log: yahs_log,
                ragtag_scaffolds: rt_scaf,
                ragtag_agp: rt_agp,
                ragtag_stats: rt_stats,
                ragtag_paf: rt_paf,
                ragtag_delta: rt_delta,
                ragtag_confidence: rt_conf,
                pretext: pretext,
                pretext_png: png,
                cool: cool,
                hic: hic,
            ]]]
        }

    emit:
    scaffolded_assemblies = ch_assemblies_scaffolded
    scaffolding_output    = ch_scaffolding_output
}
