include { FASTK_FASTK        } from "../../../modules/nf-core/fastk/fastk"
include { MERQURYFK_HAPMAKER } from "../../../modules/nf-core/merquryfk/hapmaker"
include { YAK_COUNT          } from "../../../modules/sanger-tol/yak/count"

workflow BUILD_KMER_DATABASES {
    take:
    ch_specs        // spec
    ch_data         // [meta, [reads], [fk_hist, [fk_ktabs]] ]
    val_kmer_size   // integer: kmer size

    main:
    //
    // Logic: separate input with FastK databases from input without
    //
    ch_fastk_status = ch_data
        .branch { meta, reads, fastk ->
            skip_fastk: fastk || !(meta.platform in ["pacbio_hifi", "oxford_nanopore"])
            build_fastk: true
                return [ meta, reads ]
        }

    //
    // Module: Generate FastK databases for all read sets without one
    //
    FASTK_FASTK(ch_fastk_status.build_fastk)

    //
    // Logic: Mix the FastK outputs back with the existing data
    //
    ch_data_with_fastk = ch_fastk_status.build_fastk
        .combine(FASTK_FASTK.out.hist, by: 0)
        .combine(FASTK_FASTK.out.ktab, by: 0)
        .map { meta, reads, hist, ktab ->
            def meta_new = meta + [kmer_size: val_kmer_size]
            [ meta_new, reads, [hist, ktab] ]
        }
        .mix(ch_fastk_status.skip_fastk)

    //
    // Logic: Fiter the input datasets to get those that are being used as maternal
    // or paternal datasets, and stage them for YAK building.
    //
    ch_data_for_yak = ch_specs
        .filter { spec ->
            ["maternal_dataset", "maternal_platform", "paternal_dataset", "paternal_platform"]
            .every { dataset -> spec[dataset] }
        }
        .flatMap { spec ->
            spec = spec[0]
            return [
                [
                    dataset: spec.maternal_dataset,
                    platform: spec.maternal_platform,
                ], [
                    dataset: spec.paternal_dataset,
                    platform: spec.paternal_platform
                ]
            ]
        }
        .combine(ch_data_with_fastk)
        .filter { platform_key, data_meta, _reads, _fastk ->
            data_meta.id == platform_key.dataset && data_meta.platform == platform_key.platform
        }
        .map { _platform_key, data_meta, reads, _fastk -> [data_meta, reads] }

    //
    // Module: Generate YAK kmer databases for maternal and paternal read sets
    //         for trio assembly with hifiasm
    //
    YAK_COUNT(ch_data_for_yak)

    ch_yakdbs = YAK_COUNT.out.yak
        .map { meta, yak ->
            def meta_new = meta + [kmer_size: val_kmer_size]

            return[ meta_new, yak ]
        }

    //
    // Module: Generate trio fastk databases for maternal and paternal read sets
    //         for QC with Merquryfk
    //
    ch_hapmaker_inputs = ch_specs
        .filter { spec ->
            ["maternal_dataset", "maternal_platform", "paternal_dataset", "paternal_platform"]
            .every { dataset -> spec[dataset] }
        }
        // This combines all the datasets into a list of datasets that we can map through
        .combine(ch_data_with_fastk.map { data -> [data] }.collect())
        .map { spec, datasets ->
            def out_meta = spec.subMap([
                "long_read_dataset",
                "long_read_platform",
                "maternal_dataset",
                "maternal_platform",
                "paternal_dataset",
                "paternal_platform"
            ])

            // Extract the required FastK ktabs
            def mat = datasets.find { meta, _reads, _fastk ->
                meta.id == spec.maternal_dataset && meta.platform == spec.maternal_platform
            }.get(2).get(1)
            def pat = datasets.find { meta, _reads, _fastk ->
                meta.id == spec.paternal_dataset && meta.platform == spec.paternal_platform
            }.get(2).get(1)
            def child = datasets.find { meta, _reads, _fastk ->
                meta.id == spec.long_read_dataset && meta.platform == spec.long_read_platform
            }.get(2).get(1)

            [out_meta, mat, pat, child]
        }
        .unique()
        .multiMap { meta, mat, pat, child ->
            mat: [meta, mat]
            pat: [meta, pat]
            child: [meta, child]
        }

    MERQURYFK_HAPMAKER(
        ch_hapmaker_inputs.mat,
        ch_hapmaker_inputs.pat,
        ch_hapmaker_inputs.child
    )

    ch_merqury_haptabs = MERQURYFK_HAPMAKER.out.mat_hap_ktab
        .combine(MERQURYFK_HAPMAKER.out.pat_hap_ktab, by: 0)
        .map { meta, mat, pat ->
            return meta + [mat_haptab: mat, pat_haptab: pat]
        }
        .collect()
        .map { map -> [map] }

    //
    // Logic: now that all data is pre-processed, move everything into the meta map so that each
    // data channel entry is a single map object
    //
    ch_out_data = ch_data_with_fastk
        .join(ch_yakdbs, remainder: true)
        .map { meta, reads_list, fastk, yakdb ->
            def out_meta = meta + [
                reads: reads_list.sort { f -> f.getName() },
                fk_hist: fastk[0] ?: [],
                fk_ktab: fastk[1] ?: [],
                yak: yakdb,
                haptab: []
            ]

            return out_meta
        }

    emit:
    data = ch_out_data
    merqury_trio_haptabs = ch_merqury_haptabs
}
