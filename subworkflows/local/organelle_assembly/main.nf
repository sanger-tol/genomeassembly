include { CAT_CAT as CONCATENATE_ASSEMBLIES              } from '../../../modules/nf-core/cat/cat/main'
include { MITOHIFI_FINDMITOREFERENCE                     } from '../../../modules/nf-core/mitohifi/findmitoreference/main'
include { MITOHIFI_MITOHIFI as MITOHIFI_MITOHIFI_READS   } from '../../../modules/nf-core/mitohifi/mitohifi/main'
include { MITOHIFI_MITOHIFI as MITOHIFI_MITOHIFI_CONTIGS } from '../../../modules/nf-core/mitohifi/mitohifi/main'
include { OATK                                           } from '../../../modules/nf-core/oatk/main'


workflow ORGANELLE_ASSEMBLY {
    take:
    ch_assemblies   // channel: [ val(meta), asm1, asm2  ]
    ch_long_reads   // channel: [ val(meta), long_reads  ]
    val_species     // channel: [ val(meta), val(species) ]
    val_mito_hmm    // list: [ hmm_files ]
    val_plastid_hmm // list: [ hmm_files ]

    main:
    ch_versions         = Channel.empty()
    ch_mitohifi_reads   = Channel.empty()
    ch_mitohifi_contigs = Channel.empty()

    //
    // Logic: Mitohifi does not support Conda
    //
    if(workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() == 0) {
        //
        // Module: Download reference organelle assembly
        //
        MITOHIFI_FINDMITOREFERENCE(val_species)
        ch_versions = ch_versions.mix(MITOHIFI_FINDMITOREFERENCE.out.versions)

        //
        // Module: Assemble mitogenome from reads using MitoHiFi
        //
        ch_mitohifi_reads_input = ch_long_reads
            | combine(MITOHIFI_FINDMITOREFERENCE.out.fasta)
            | combine(MITOHIFI_FINDMITOREFERENCE.out.gb)
            | multiMap { meta, reads, meta_fasta, fasta, meta_gb, gb ->
                reads: [ meta, reads ]
                fasta: [ meta, fasta ]
                gb:    [ meta, gb ]
                method: "reads"
                code: meta.code
            }

        MITOHIFI_MITOHIFI_READS(
            ch_mitohifi_reads_input.reads,
            ch_mitohifi_reads_input.fasta,
            ch_mitohifi_reads_input.gb,
            ch_mitohifi_reads_input.method,
            ch_mitohifi_reads_input.code
        )
        ch_versions = ch_versions.mix(MITOHIFI_MITOHIFI_READS.out.versions)

        //
        // Module: Concatenate assembly pairs for Mitohifi
        //
        ch_concat_input = ch_assemblies
            | map { meta, asm1, asm2 -> [ meta, [asm1, asm2] ] }

        CONCATENATE_ASSEMBLIES(ch_concat_input)
        ch_versions = ch_versions.mix(CONCATENATE_ASSEMBLIES.out.versions)

        //
        // Module: Identify organelle from assembled contigs with Mitohifi
        //
        ch_mitohifi_contigs_input = ch_assemblies
            | combine(MITOHIFI_FINDMITOREFERENCE.out.fasta)
            | combine(MITOHIFI_FINDMITOREFERENCE.out.gb)
            | multiMap { meta, asm, meta_fasta, fasta, meta_gb, gb ->
                asm: [ meta, reads ]
                fasta: [ meta, fasta ]
                gb:    [ meta, gb ]
                method: "contigs"
                code: meta.code
            }
        MITOHIFI_MITOHIFI_CONTIGS(
            ch_mitohifi_contigs_input.asm,
            ch_mitohifi_contigs_input.fasta,
            ch_mitohifi_contigs_input.gb,
            ch_mitohifi_contigs_input.method,
            ch_mitohifi_contigs_input.code
        )
        ch_versions = ch_versions.mix(MITOHIFI_MITOHIFI_CONTIGS.out.versions)

        //
        // Logic: Prepare all outputs from Mitohifi for emission
        //        Do it this way as we will move to a channel publishing structure in future
        //
        ch_mitohifi_reads = Channel.empty()
            | mix(
                MITOHIFI_MITOHIFI_READS.out.fasta,
                MITOHIFI_MITOHIFI_READS.out.stats,
                MITOHIFI_MITOHIFI_READS.out.gb,
                MITOHIFI_MITOHIFI_READS.out.gff,
                MITOHIFI_MITOHIFI_READS.out.all_potential_contigs,
                MITOHIFI_MITOHIFI_READS.out.contigs_annotations,
                MITOHIFI_MITOHIFI_READS.out.contigs_circularization,
                MITOHIFI_MITOHIFI_READS.out.contigs_filtering,
                MITOHIFI_MITOHIFI_READS.out.coverage_mapping,
                MITOHIFI_MITOHIFI_READS.out.coverage_plot,
                MITOHIFI_MITOHIFI_READS.out.final_mitogenome_annotation,
                MITOHIFI_MITOHIFI_READS.out.final_mitogenome_choice,
                MITOHIFI_MITOHIFI_READS.out.final_mitogenome_coverage,
                MITOHIFI_MITOHIFI_READS.out.potential_contigs,
                MITOHIFI_MITOHIFI_READS.out.reads_mapping_and_assembly,
                MITOHIFI_MITOHIFI_READS.out.shared_genes,
            )

        ch_mitohifi_contigs = Channel.empty()
            | mix(
                MITOHIFI_MITOHIFI_READS.out.fasta,
                MITOHIFI_MITOHIFI_READS.out.stats,
                MITOHIFI_MITOHIFI_READS.out.gb,
                MITOHIFI_MITOHIFI_READS.out.gff,
                MITOHIFI_MITOHIFI_READS.out.all_potential_contigs,
                MITOHIFI_MITOHIFI_READS.out.contigs_annotations,
                MITOHIFI_MITOHIFI_READS.out.contigs_circularization,
                MITOHIFI_MITOHIFI_READS.out.contigs_filtering,
                MITOHIFI_MITOHIFI_READS.out.coverage_mapping,
                MITOHIFI_MITOHIFI_READS.out.coverage_plot,
                MITOHIFI_MITOHIFI_READS.out.final_mitogenome_annotation,
                MITOHIFI_MITOHIFI_READS.out.final_mitogenome_choice,
                MITOHIFI_MITOHIFI_READS.out.final_mitogenome_coverage,
                MITOHIFI_MITOHIFI_READS.out.potential_contigs,
                MITOHIFI_MITOHIFI_READS.out.reads_mapping_and_assembly,
                MITOHIFI_MITOHIFI_READS.out.shared_genes,
            )
    }

    //
    // Module: Assemble mito from reads using Oatk
    //
    OATK(
        ch_long_reads,
        val_mito_hmm,
        val_plastid_hmm
    )
    ch_versions = ch_versions.mix(OATK.out.versions)

    ch_oatk_out_mito = Channel.empty()
        | mix(
            OATK.out.mito_fasta,
            OATK.out.pltd_fasta,
            OATK.out.mito_bed,
            OATK.out.pltd_bed,
            OATK.out.mito_gfa,
            OATK.out.pltd_gfa,
            OATK.out.annot_mito_txt,
            OATK.out.annot_pltd_txt,
            OATK.out.clean_gfa,
            OATK.out.final_gfa,
            OATK.out.initial_gfa,
            OATK.out.multiplex_gfa,
            OATK.out.unzip_gfa,
        )

    emit:
    mito_mitohifi_reads   = ch_mitohifi_reads
    mito_mitohifi_contigs = ch_mitohifi_contigs
    organelles_oatk       = ch_oatk_out
    versions              = ch_versions
}
