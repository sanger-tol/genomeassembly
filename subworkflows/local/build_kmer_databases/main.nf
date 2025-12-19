include { FASTK_FASTK        } from "../../../modules/nf-core/fastk/fastk/main"
include { YAK_COUNT          } from "../../../modules/nf-core/yak/count/main"
include { MERQURYFK_HAPMAKER } from "../../../modules/nf-core/merquryfk/hapmaker"

workflow BUILD_KMER_DATABASES {
    take:
    ch_long_reads   // [meta, [reads], fk_hist, [fk_ktab]]
    ch_trio         // [trio_meta, [mat_meta, [reads], fk_hist, [fk_ktab]], [pat_meta, [reads], fk_hist, [fk_ktab]], [child_meta, [reads], fk_hist, [fk_ktab]], ]
    val_kmer_size   // integer: kmer size

    main:
    ch_versions = channel.empty()

    //
    // Module: Generate FastK databases for all read sets without one
    //
    ch_fastk_status = ch_long_reads.mix(ch_trio.transpose())
        .branch { meta, reads, hist, _ktab ->
            has_fastk: hist
            no_fastk: true
                return [ meta, reads ]
        }

    FASTK_FASTK(ch_fastk_status.no_fastk)
    ch_versions = ch_versions.mix(FASTK_FASTK.out.versions)

    ch_fastk = ch_fastk_status.no_fastk
        .combine(FASTK_FASTK.out.hist, by: 0)
        .combine(FASTK_FASTK.out.ktab, by: 0)
        .map { meta, reads, hist, ktab ->
            def meta_new = meta + [kmer_size: val_kmer_size]
            [meta_new, reads, hist, ktab]
        }
        .mix(ch_fastk_status.has_fastk)

    ch_long_reads = ch_fastk.filter { meta, _reads, _hist, _ktab -> meta.read_type == "long" }

    //
    // Module: Generate trio databases for maternal and paternal read sets
    //         for trio assembly with hifiasm
    //
    ch_yak_input = ch_trio.transpose()
        .filter { meta, _reads, _hist, _ktab -> meta.trio in ["mat", "pat"]}
        .map { meta, reads, _hist, _ktab -> [meta, reads] }

    YAK_COUNT(ch_yak_input)
    ch_versions = ch_versions.mix(YAK_COUNT.out.versions)

    ch_trio_yak_dbs = YAK_COUNT.out.yak
        .map { meta, yak -> [ meta.trio_id, meta, yak ] }
        .groupTuple(by: 0, size: 2)
        .map { trio_id, _meta, yaks ->
            def pat = yaks.find { yak -> yak.name =~ /pat.yak$/ }
            def mat = yaks.find { yak -> yak.name =~ /mat.yak$/ }
            [ [id: trio_id], pat, mat ]
        }

    //
    // Module: Generate trio fastk databases for maternal and paternal read sets
    //         for QC with Merquryfk
    //
    ch_trio_ktabs_split = ch_fastk
        .filter { meta, _reads, _hist, _ktab -> meta.read_type == "trio" } 
        .branch { meta, _reads, _hist, ktab -> 
            mat: meta.trio == "mat"
                return [ meta - meta.subMap("trio"), ktab ]
            pat: meta.trio == "pat"
                return [ meta - meta.subMap("trio"), ktab ]
            child: meta.trio == "child"
                return [ meta - meta.subMap("trio"), ktab ]
        }

    ch_hapmaker_input = ch_trio_ktabs_split.mat
        .join(ch_trio_ktabs_split.pat, by: 0)
        .join(ch_trio_ktabs_split.child, by: 0)
        .multiMap { meta, mat, pat, child -> 
            mat: [meta, mat]
            pat: [meta, pat]
            child: [meta, child]
        }

    MERQURYFK_HAPMAKER(
        ch_hapmaker_input.mat,
        ch_hapmaker_input.pat,
        ch_hapmaker_input.child
    )
    ch_versions = ch_versions.mix(MERQURYFK_HAPMAKER.out.versions)

    ch_hapdbs = MERQURYFK_HAPMAKER.out.mat_hap_ktab.join(MERQURYFK_HAPMAKER.out.pat_hap_ktab)

    emit:
    long_reads     = ch_long_reads
    trio_yakdb     = ch_trio_yak_dbs
    trio_hapdbs    = ch_hapdbs
    versions       = ch_versions
}
