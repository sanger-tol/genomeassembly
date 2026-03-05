include { MITOHIFI_FINDMITOREFERENCE        } from '../../../modules/nf-core/mitohifi/findmitoreference'
include { MITOHIFI_MITOHIFI                 } from '../../../modules/nf-core/mitohifi/mitohifi'

workflow MITOHIFI_ASSEMBLY {
    take:
    ch_mitohifi_specs  // channel: spec
    ch_assemblies      // channel: [spec, hap1, hap2]

    main:
    //
    // Logic: Mitohifi does not support Conda
    //
    if(workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() == 0) {
        //
        // Logic: set up the inputs to mitohifi findMitoReference script to download
        // reference mito and plastid genomes for reference-based organelle assembly
        //
        ch_species_to_download = ch_mitohifi_specs
            .map { spec ->
                [spec.params.subMap(["mitohifi_reference_species", "organelle"]), spec.params.mitohifi_reference_species]
            }
            .unique()

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
                reads_input: spec.params.mode == "reads"
                assembly_input: spec.params.mode == "contigs"
            }

        // Stage the mitohifi assembly input. We need to match the assemblies
        // to the spec by matching the spec hash to the assembly hash
        ch_mitohifi_asm_input = ch_mitohifi_specs_split.assembly_input
            .combine(ch_assemblies)
            .combine(MITOHIFI_FINDMITOREFERENCE.out.reference)
            .filter { spec, asm_meta, _asm1, _asm2, ref_meta, _ref_fa, _ref_gb ->
                def hash_match = spec.prevID == asm_meta.id
                def species_match = spec.params.mitohifi_reference_species == ref_meta.mitohifi_reference_species
                def organelle_match = spec.params.organelle == ref_meta.organelle

                hash_match && species_match && organelle_match
            }
            .map { spec, _asm_meta, asm1, asm2, _ref_meta, ref_fa, ref_gb ->
                def params_out = spec.params + [
                    mitohifi_reference_fa: ref_fa,
                    mitohifi_reference_gb: ref_gb,
                ]
                return [spec + params_out, [asm1, asm2], ref_fa, ref_gb]
            }

        ch_mitohifi_reads_input = ch_mitohifi_specs_split.reads_input
            .combine(MITOHIFI_FINDMITOREFERENCE.out.reference)
            .filter { spec, ref_meta, _ref_fa, _ref_gb ->
                def species_match = spec.params.mitohifi_reference_species == ref_meta.mitohifi_reference_species
                def organelle_match = spec.params.organelle == ref_meta.organelle

                species_match && organelle_match
            }
            .map { spec, _ref_meta, ref_fa, ref_gb ->
                def params_out = spec.params + [
                    mitohifi_reference_fa: ref_fa,
                    mitohifi_reference_gb: ref_gb,
                ]
                return [spec + params_out, spec.data.long_read.reads, ref_fa, ref_gb]
            }


        ch_mitohifi_input = ch_mitohifi_asm_input.mix(ch_mitohifi_reads_input)
            .multiMap { spec, input, ref_fa, ref_gb ->
                def genetic_code = spec.params.organelle == "mito"
                    ? spec.params.mitohifi_mito_genetic_code
                    : spec.params.mitohifi_plastid_genetic_code

                input: [spec, input]
                reference: [spec, ref_fa, ref_gb]
                method: spec.params.mode
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
        ch_mitohifi_output = MITOHIFI_MITOHIFI.out.all_files
            .map { spec, mitohifi_files ->

                // The potential contigs directory contains symlinks which cannot be published.
                // Loop through all contents and discard these.
                def valid_publish_files = []
                mitohifi_files.each { mitohifi_file ->
                    if(mitohifi_file.isDirectory()) {
                        mitohifi_file.eachFileRecurse { contents ->
                            if(!contents.isLink() && !contents.isDirectory()) {
                                valid_publish_files << contents
                            }
                        }
                    } else if(!mitohifi_file.isLink()) {
                        valid_publish_files << mitohifi_file
                    }
                }

                return spec.subMap(["id", "stage", "data", "params", "tools"]) + [
                    output: [
                        mitohifi: [
                            mitohifi_reference_fa: spec.params.mitohifi_reference_fa,
                            mitohifi_reference_gb: spec.params.mitohifi_reference_gb,
                            mitohifi_files: valid_publish_files
                        ]
                    ]
                ]
            }
    }

    emit:
    mitohifi_assemblies = ch_mitohifi_output
}
