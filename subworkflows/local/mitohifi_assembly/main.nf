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
            .filter { spec -> spec.params?.download_mito_reference || spec.params?.download_plastid_reference }
            .map { spec ->
                [spec.params.subMap(["mitohifi_reference_species", "organelle"]), spec.params.mitohifi_reference_species]
            }
            .unique()

        //
        // Module: Download reference organelle assembly
        //
        ch_mitohifi_references = MITOHIFI_FINDMITOREFERENCE(ch_species_to_download).reference
            .map { meta, ref_fa, ref_gb -> meta + [ref_fa: ref_fa, ref_gb: ref_gb] }
            .ifEmpty([organelle: null, mitohifi_reference_species: null, ref_fa: null, ref_gb: null])

        //
        // Logic: Separate out the two types of mitohifi input to allow safe joining
        // without worrying about empty channels.
        //
        ch_mitohifi_specs_split = ch_mitohifi_specs
            .branch { spec ->
                reads_input: spec.params.mode == "reads"
                assembly_input: spec.params.mode == "contigs"
            }

        //
        // Logic: Stage the mitohifi assembly input. For contigs mode, we need to match
        // the assemblies to the spec by matching the spec hash to the assembly hash.
        //
        // For the reads, we just need to pull the reads out into a new tuple entry.
        //
        ch_mitohifi_asm_input = ch_mitohifi_specs_split.assembly_input
            .combine(ch_assemblies)
            .filter { spec, asm_meta, _asm1, _asm2 -> spec.prevID == asm_meta.id }
            .map { spec, _asm_meta, asm1, asm2 -> [spec, [asm1, asm2]] }

        ch_mitohifi_reads_input = ch_mitohifi_specs_split.reads_input
            .map { spec -> [spec, spec.data.long_read.reads] }

        //
        // Logic: Set up the reference files for mitohifi. Combine with a collected list of
        // all downloaded references, and if the ref_fa and ref_gb parts are not set, then
        // find the requisite species/organelle reference.
        //
        // Also drops the original keys and replaces them with a more generic keyset.
        //
        ch_mitohifi_input = ch_mitohifi_asm_input.mix(ch_mitohifi_reads_input)
            .combine(ch_mitohifi_references.collect().map { refs -> [refs] })
            .map { spec, input, references ->
                def organelle = spec.params.organelle

                def required_reference = references.find { ref ->
                    ref.organelle == organelle && spec.params.mitohifi_reference_species == ref.mitohifi_reference_species
                }

                def ref_fa = spec.params["mitohifi_${organelle}_reference_fa"] ?: required_reference.ref_fa
                def ref_gb = spec.params["mitohifi_${organelle}_reference_gb"] ?: required_reference.ref_gb

                // remove old keys
                def keysToRemove = [
                    "mitohifi_mito_reference_fa",
                    "mitohifi_mito_reference_gb",
                    "mitohifi_plastid_reference_fa",
                    "mitohifi_plastid_reference_gb"
                ]

                def params_out = (spec.params - spec.params.subMap(keysToRemove)) + [
                    mitohifi_reference_fa: ref_fa,
                    mitohifi_reference_gb: ref_gb,
                ]

                return [spec + [params: params_out], input, ref_fa, ref_gb]
            }.view()
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
                def downloaded_reference = spec.params.download_mito_reference || spec.params.download_plastid_reference

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
                            mitohifi_reference_fa: downloaded_reference ? spec.params.mitohifi_reference_fa : null,
                            mitohifi_reference_gb: downloaded_reference ? spec.params.mitohifi_reference_gb : null,
                            mitohifi_files: valid_publish_files
                        ]
                    ]
                ]
            }
    }

    emit:
    mitohifi_assemblies = ch_mitohifi_output
}
