include { OATK_ASSEMBLY     } from '../../../subworkflows/local/oatk_assembly'
include { MITOHIFI_ASSEMBLY } from '../../../subworkflows/local/mitohifi_assembly'

workflow ORGANELLE_ASSEMBLY {
    take:
    ch_organelle_specs

    main:
    ch_versions = channel.empty()

    ch_oatk_specs = ch_organelle_specs.filter { spec -> spec.assembler == "oatk" }.map { spec -> spec.oatk }.unique()
    ch_mitohifi_specs = ch_organelle_specs.filter { spec -> spec.assembler == "mitohifi" }.map { spec -> spec.mitohifi }.unique()

    //
    // Subworkflow: assemble organelles with oatk
    //
    OATK_ASSEMBLY(
        ch_oatk_specs
    )

    //
    // Subworkflow: assemble organelles with mitohifi
    //
    MITOHIFI_ASSEMBLY(
        ch_mitohifi_specs,
        channel.empty(),
    )

    emit:
    ch_oatk_output = OATK_ASSEMBLY.out.oatk_output
    ch_mitohifi_output = MITOHIFI_ASSEMBLY.out.mitohifi_assemblies
}
