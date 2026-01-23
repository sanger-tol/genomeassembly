include { MITOHIFI_FINDMITOREFERENCE        } from '../../../modules/nf-core/mitohifi/findmitoreference/main'
include { MITOHIFI_MITOHIFI                 } from '../../../modules/nf-core/mitohifi/mitohifi/main'

workflow MITOHIFI {
    take:
    ch_mitohifi_specs  // channel: spec
    ch_assemblies      // channel: [spec, hap1, hap2]
    ch_data            // channel: data_map

    main:
    //
    // Logic: set up the inputs to mitohifi findMitoReference script to download
    // reference mito and plastid genomes for reference-based organelle assembly
    //
    ch_species_to_download = ch_mitohifi_specs
        .map { spec ->
            [spec.subMap(["mitohifi_reference_species", "organelle"]), spec.mitohifi_reference_species]
        }
        .unique()

    //
    // Logic: Mitohifi does not support Conda
    //
    if(workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() == 0) {
        //
        // Module: Download reference organelle assembly
        //
        MITOHIFI_FINDMITOREFERENCE(ch_species_to_download)

        //
        // Logic: Separate out the two types of mitohifi input to allow safe joining
        // without worrying about empty channels.
        //
        ch_mitohifi_specs_split = ch_mitohifi_specs
            .branch { spec ->
                reads_input: spec.mode == "reads"
                assembly_input: spec.mode == "contigs"
            }

        // Stage the mitohifi assembly input. We need to match the assemblies
        // to the spec by matching the spec hash to the assembly hash
        ch_mitohifi_asm_input = ch_mitohifi_specs_split.assembly_input
            .combine(ch_assemblies)
            .combine(MITOHIFI_FINDMITOREFERENCE.out.reference)
            .filter { spec, asm_meta, _asm1, _asm2, ref_meta, _ref_fa, _ref_gb ->
                def hash_match = spec.prevHash == asm_meta.hash
                def species_match = spec.mitohifi_reference_species == ref_meta.mitohifi_reference_species
                def organelle_match = spec.organelle == ref_meta.organelle

                hash_match && species_match && organelle_match
            }
            .map { spec, _asm_meta, asm1, asm2, _ref_meta, ref_fa, ref_gb ->
                return [spec, [asm1, asm2], ref_fa, ref_gb]
            }

        ch_mitohifi_reads_input = ch_mitohifi_specs_split.reads_input
            .combine(ch_data)
            .combine(MITOHIFI_FINDMITOREFERENCE.out.reference)
            .filter { spec, datasets, ref_meta, _ref_fa, _ref_gb ->
                def species_match = spec.mitohifi_reference_species == ref_meta.mitohifi_reference_species
                def organelle_match = spec.organelle == ref_meta.organelle

                species_match && organelle_match
            }
            .map { spec, datasets, _ref_meta, ref_fa, ref_gb ->
                def long_reads = datasets.find { data ->
                    data.id == spec.long_read_dataset && data.platform == spec.long_read_platform
                }

                return [spec, long_reads.reads, ref_fa, ref_gb]
            }


        ch_mitohifi_input = ch_mitohifi_asm_input.mix(ch_mitohifi_reads_input)
            .multiMap { spec, input, ref_fa, ref_gb ->
                def genetic_code = spec.organelle == "mito" ? spec.mitohifi_genetic_code : spec.mitohifi_plastid_genetic_code

                input: [ spec, input ]
                reference: [ spec, ref_fa, ref_gb ]
                method: spec.mode
                code: genetic_code
            }

        //
        // Module: Assemble or extract mitogenome using MitoHiFi
        //
        MITOHIFI_MITOHIFI(
            ch_mitohifi_input.input,
            ch_mitohifi_input.reference,
            ch_mitohifi_input.method,
            ch_mitohifi_input.code
        )

        //
        // Logic: Prepare all outputs from Mitohifi for emission
        //        Do it this way as we will move to a channel publishing structure in future
        //
        ch_mitohifi_output = MITOHIFI_MITOHIFI.out.fasta
                .join(MITOHIFI_MITOHIFI.out.stats, by: 0)
                .join(MITOHIFI_MITOHIFI.out.gb, by: 0, remainder: true)
                .join(MITOHIFI_MITOHIFI.out.gff, by: 0, remainder: true)
                .join(MITOHIFI_MITOHIFI.out.all_potential_contigs, by: 0)
                .join(MITOHIFI_MITOHIFI.out.contigs_annotations, by: 0)
                .join(MITOHIFI_MITOHIFI.out.contigs_circularization, by: 0)
                .join(MITOHIFI_MITOHIFI.out.contigs_filtering, by: 0)
                .join(MITOHIFI_MITOHIFI.out.coverage_mapping, by: 0, remainder: true)
                .join(MITOHIFI_MITOHIFI.out.coverage_plot, by: 0, remainder: true)
                .join(MITOHIFI_MITOHIFI.out.final_mitogenome_annotation, by: 0)
                .join(MITOHIFI_MITOHIFI.out.final_mitogenome_choice, by: 0)
                .join(MITOHIFI_MITOHIFI.out.final_mitogenome_coverage, by: 0)
                .join(MITOHIFI_MITOHIFI.out.potential_contigs, by: 0)
                .join(MITOHIFI_MITOHIFI.out.reads_mapping_and_assembly, by: 0, remainder: true)
                .join(MITOHIFI_MITOHIFI.out.shared_genes, by: 0)
                .join(MITOHIFI_MITOHIFI.out.log, by: 0)
            .map { spec, fasta, gb, gff, apc, annot, circ, filt, map, cov_plot, final_annot,
                final_choice, final_coverage, pot_contigs, rma, shared, log ->
                return [
                    id: spec.hash,
                    organelle: spec.organelle,
                    fasta: fasta,
                    stats: stats,
                    gb: gb,
                    gff: gff,
                    all_potential_contigs: apc,
                    contigs_annotations: annot,
                    contigs_circularisation: circ,
                    contigs_filtering: filt,
                    coverage_mapping: map,
                    coverage_plot: cov_plot,
                    final_mitogenome_annotation: final_annot,
                    final_mitogenome_choice: final_choice,
                    final_mitogenome_coverage: final_coverage,
                    potential_contigs: pot_contigs,
                    reads_mapping_and_assembly: rma,
                    shared_genes: shared_genes,
                    log: log
                ]
            }
    }

    emit:
    assemblies = ch_mitohifi_output
}
