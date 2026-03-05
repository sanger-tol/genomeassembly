include { CRAM_MAP_ILLUMINA_HIC           } from '../../../subworkflows/sanger-tol/cram_map_illumina_hic'
include { BAM_STATS_SAMTOOLS              } from '../../../subworkflows/nf-core/bam_stats_samtools'
include { FASTA_BAM_SCAFFOLDING_YAHS      } from '../../../subworkflows/sanger-tol/fasta_bam_scaffolding_yahs'

include { TABIX_BGZIP as BGZIP_SCAFFOLDED } from '../../../modules/nf-core/tabix/bgzip'

workflow SCAFFOLDING {
    take:
    ch_scaffolding_specs            // spec
    ch_assemblies                   // [meta, hap1, hap2]
    val_hic_aligner                 // "bwamem2" or "minimap2"
    val_hic_mapping_cram_chunk_size // int > 1
    val_cool_bin                    // int > 1

    main:
    //
    // Logic: join all the assemblies with the scaffolding specifications and
    // data, filter for those assemblies which are to be scaffolded, then
    // map out the data for the purging subworkflow.
    //
    ch_hic_mapping_inputs = ch_assemblies
        .combine(ch_scaffolding_specs)
        .filter { asm_meta, _asm1, _asm2, spec -> asm_meta.id == spec.prevID }
        .multiMap { _asm_meta, asm1, asm2, spec ->
            def spec_hap1 = spec + [_hap: "hap1"]
            def spec_hap2 = spec + [_hap: "hap2"]
            hap1: [ spec_hap1, asm1 ]
            hap2: [ spec_hap2, asm2 ]
            hic_reads: [ [spec_hap1, spec_hap2], spec.data.hic.reads ]
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
    // Subworkflow: Calculate stats for Hi-C mapping
    //
    ch_hic_mapping_stats_input = CRAM_MAP_ILLUMINA_HIC.out.bam
        .combine(CRAM_MAP_ILLUMINA_HIC.out.bam_index.filter { _meta, idx -> idx.getExtension() == "csi" }, by: 0)
        .combine(ch_hic_mapping_inputs.hap1.mix(ch_hic_mapping_inputs.hap2), by: 0)
        .multiMap { meta, bam, bai, asm ->
            bam: [ meta, bam, bai ]
            asm: [ meta, asm ]
        }

    BAM_STATS_SAMTOOLS(
        ch_hic_mapping_stats_input.bam,
        ch_hic_mapping_stats_input.asm
    )

    //
    // Subworkflow: scaffold assemblies using yahs and create contact maps
    //
    FASTA_BAM_SCAFFOLDING_YAHS(
        ch_hic_mapping_inputs.hap1.mix(ch_hic_mapping_inputs.hap2),
        CRAM_MAP_ILLUMINA_HIC.out.bam,
        true,
        true,
        true,
        val_cool_bin
    )

    //
    // Module: bgzip all scaffolded assembly fasta
    //
    BGZIP_SCAFFOLDED(FASTA_BAM_SCAFFOLDING_YAHS.out.scaffolds_fasta)

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
        .join(CRAM_MAP_ILLUMINA_HIC.out.bam, by: 0)
        .join(CRAM_MAP_ILLUMINA_HIC.out.bam_index.filter { _meta, idx -> idx.getExtension() == "csi" }, by: 0)
        .join(BAM_STATS_SAMTOOLS.out.stats, by: 0)
        .join(BAM_STATS_SAMTOOLS.out.flagstat, by: 0)
        .join(BAM_STATS_SAMTOOLS.out.idxstats, by: 0)
        .join(FASTA_BAM_SCAFFOLDING_YAHS.out.scaffolds_agp, by: 0)
        .join(FASTA_BAM_SCAFFOLDING_YAHS.out.yahs_bin, by: 0)
        .join(FASTA_BAM_SCAFFOLDING_YAHS.out.yahs_inital, by: 0, remainder: true)
        .join(FASTA_BAM_SCAFFOLDING_YAHS.out.yahs_intermediate, by: 0, remainder: true)
        .join(FASTA_BAM_SCAFFOLDING_YAHS.out.yahs_log, by: 0)
        .join(FASTA_BAM_SCAFFOLDING_YAHS.out.pretext, by: 0, remainder: true)
        .join(FASTA_BAM_SCAFFOLDING_YAHS.out.pretext_png, by: 0, remainder: true)
        .join(FASTA_BAM_SCAFFOLDING_YAHS.out.cool, by: 0, remainder: true)
        .join(FASTA_BAM_SCAFFOLDING_YAHS.out.hic, by: 0, remainder: true)
        .map { spec, fasta, bam, bai, stats, flagstats, idxstats, agp, bin, initial, intermed, log, pretext, png, cool, hic ->
            return spec.subMap(["id", "stage", "data", "params", "tools"]) + [
                hap: spec._hap,
                output: [
                    scaffolding: [
                        fasta: fasta,
                        bam: bam,
                        bai: bai,
                        stats: stats,
                        flagstats: flagstats,
                        idxstats: idxstats,
                        yahs_agp: agp,
                        yahs_bin: bin,
                        yahs_initial: initial,
                        yahs_intermeriate: intermed,
                        yahs_log: log,
                        pretext: pretext,
                        pretext_png: png,
                        cool: cool,
                        hic: hic
                    ]
                ]
            ]
        }

    emit:
    scaffolded_assemblies = ch_assemblies_scaffolded
    scaffolding_output = ch_scaffolding_output
}
