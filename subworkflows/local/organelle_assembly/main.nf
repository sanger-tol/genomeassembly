include { OATK_ASSEMBLY     } from '../../../subworkflows/local/oatk_assembly'
include { MITOHIFI_ASSEMBLY } from '../../../subworkflows/local/mitohifi_assembly'

include { setupStage        } from '../../../functions/local/assembly_stages'

include { GENERATE_SPECIFICATION_INDEX as INDEX_SPEC } from '../../../modules/local/generate_specification_index'

workflow ORGANELLE_ASSEMBLY {
    take:
    ch_specs

    main:
    //
    // Logic: Set up and deduplicate input stages
    //
    ch_stages = ch_specs.flatMap { spec ->
        spec.stages.collect { stageName, _stageData -> setupStage(spec, stageName) }
    }.unique()

    //
    // Subworkflow: assemble organelles with Oatk
    //
    OATK_ASSEMBLY(
        ch_stages.filter { stage -> stage.stage == "oatk" }
    )

    //
    // Subworkflow: assemble organelles with Mitohifi
    //
    MITOHIFI_ASSEMBLY(
        ch_stages.filter { stage -> stage.stage == "mitohifi" },
        channel.empty(),
    )

    //
    // Module: Generate index files for each stage and the overall stage
    //
    ch_versions_to_index = channel.topic("versions")
        .map { task, tool, version -> [task.split(':').last(), tool, version] }
        .unique()
        .groupTuple(by: 0)
        .map { task, tools, versions ->
            [(task): [tools, versions].transpose().collectEntries()]
        }
        .reduce { a, b -> a + b }

    INDEX_SPEC(
        ch_specs.map { spec -> spec.subMap(["name", "assembler", "data", "params", "tools"]) },
        ch_versions_to_index
    )

    //
    // Logic: Re-join stages to their input specifications for publishing
    //
    ch_output_specs = ch_specs.map { spec -> spec.subMap(["name", "hashes", "assembler"]) }.unique()

    ch_oatk_output = ch_output_specs
        .combine(OATK_ASSEMBLY.out.oatk_output)
        .filter { spec, oatk -> oatk.id in spec.hashes.values() }
        .map { spec, oatk -> spec + oatk }

    ch_mitohifi_output = ch_output_specs
        .combine(MITOHIFI_ASSEMBLY.out.mitohifi_assemblies)
        .filter { spec, mitohifi -> mitohifi.id in spec.hashes.values() }
        .map { spec, mitohifi -> spec + mitohifi }

    emit:
    oatk           = ch_oatk_output
    reads_mitohifi = ch_mitohifi_output
    spec_indexes   = INDEX_SPEC.out.spec
}
