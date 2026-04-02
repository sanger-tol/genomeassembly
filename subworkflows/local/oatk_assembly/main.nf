include { OATK          } from '../../../modules/nf-core/oatk'
include { BANDAGE_IMAGE } from '../../../modules/nf-core/bandage/image'

workflow OATK_ASSEMBLY {
    take:
    ch_oatk_specs

    main:
    //
    // Logic: join all the the oatk specifications with the data
    //
    ch_oatk_input = ch_oatk_specs
        .multiMap { spec ->
            reads: [spec, spec.data.long_read.reads]
            mito_hmm: spec.params.oatk_mito_hmm ? [spec, spec.params.oatk_mito_hmm] : [[], []]
            plastid_hmm: spec.params.oatk_plastid_hmm ? [spec, spec.params.oatk_plastid_hmm] : [[], []]
        }

    //
    // Module: assemble organelles with oatk
    //
    OATK(
        ch_oatk_input.reads,
        ch_oatk_input.mito_hmm,
        ch_oatk_input.plastid_hmm
    )

    //
    // Module: generate plots of the assembly graphs with Bandage.
    // We filter out empty files here as they cause a Bandage failure.
    //
    BANDAGE_IMAGE(
        OATK.out.mito_gfa.mix(OATK.out.pltd_gfa, OATK.out.final_gfa)
            .filter { _meta, gfa -> gfa.size() > 0 }
    )

    //
    // Logic: combine all oatk outputs into a single map for ease of publishing
    //
    ch_oatk_outputs = channel.empty()
        .mix(
            OATK.out.mito_fasta,
            OATK.out.mito_bed,
            OATK.out.mito_gfa,
            OATK.out.annot_mito_txt,
            OATK.out.pltd_fasta,
            OATK.out.pltd_bed,
            OATK.out.pltd_gfa,
            OATK.out.annot_pltd_txt,
            OATK.out.initial_gfa,
            OATK.out.final_gfa,
            OATK.out.log
        )
        .groupTuple(by: 0)
        .join(
            BANDAGE_IMAGE.out.png.groupTuple(by: 0)
        )
        .map { spec, oatk_outputs, images ->
            return spec.subMap(["id", "stage", "data", "params", "tools"]) + [
                output: [
                    oatk: [
                        oatk: oatk_outputs,
                        bandage_plots: images
                    ]
                ]
            ]
        }

    emit:
    oatk_output = ch_oatk_outputs
}
