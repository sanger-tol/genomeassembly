include { POLISHING_10X } from '../../../subworkflows/local/polishing_10X'

include { CAT_CAT as CONCATENATE_ASSEMBLIES } from '../../../modules/nf-core/cat/cat'
include { SEQKIT_GREP as SEPARATE_HAPLOTYPES } from '../../../modules/nf-core/seqkit/grep/main'
include { TABIX_BGZIP as BGZIP_POLISHED } from '../../../modules/nf-core/tabix/bgzip/main'

workflow POLISHING {
    take:
    ch_polishing_specs // channel: spec
    ch_assemblies      // channel: [spec, hap1, hap2]
    ch_data            // channel: data_map

    main:
    ch_versions = channel.empty()

    //
    // Logic: join all the assemblies with the purging specifications and
    // data, filter for those assemblies which are to be purged, then
    // map out the data for the purging subworkflow.
    //
    ch_polishing_input = ch_assemblies
        .combine(ch_polishing_specs)
        .combine(ch_data)
        .filter { meta, asm1, asm2, spec, datasets ->
            meta.hash == spec.prevHash
        }
        .multiMap { meta, asm1, asm2, spec, datasets ->
            def i10x_reads = datasets.find { data ->
                data.id == spec.illumina_10x_dataset && data.platform == "illumina_10x"
            }

            assemblies: [ spec, asm1, asm2 ]
            illumina_10x_reads: [ spec, long_reads.reads ]
        }

    //
    // Module: Concatenate hap1/hap2 together for polishing
    //
    CONCATENATE_ASSEMBLIES(ch_polishing_input.assemblies)
    ch_versions = ch_versions.mix(CONCATENATE_ASSEMBLIES.out.versions)

    //
    // Subworkflow: run polishing pipeline
    //
    POLISHING_10X(
        CONCATENATE_ASSEMBLIES.out.file_out,
        ch_polishing_input.illumina_10x_reads
    )
    ch_versions = ch_versions.mix(POLISHING_10X.out.versions)

    //
    // Module: Separate back out primary/alt/hap1/hap2 contigs
    //
    ch_polished_assemblies_to_separate = POLISHING_10X.out.polished_assemblies
        .flatMap { meta, asm, hap ->
            return (1..2).collect { count ->
                [meta + [_hap: "$hap{count}"], asm]
            }
        }

    //
    // Module: split the haplotypes back out using seqkit grep + regex
    //
    SEPARATE_HAPLOTYPES(
        ch_polished_assemblies_to_separate,
        []
    )
    ch_versions = ch_versions.mix(SEPARATE_HAPLOTYPES.out.versions)

    //
    // Logic: Take the split assemblies and re-orgnaise into a [meta, hap1, hap2] format
    //
    ch_assemblies_polished_split = SEPARATE_HAPLOTYPES.out.filter
        .map { meta, asm ->
            def meta_new = meta - meta.subMap("_hap")
            [ meta_new, meta._hap, asm ]
        }
        .groupTuple(by: 0, size: 2, remainder: true)
        .map { meta_new, hap_labels, asms ->
            [ meta_new, asms[hap_labels.indexOf("hap1")], asms[hap_labels.indexOf("hap2")] ]
        }

    //
    // Module: bgzip all polished assembly fasta
    //
    ch_bgzip_input = ch_assemblies_polished_split
        .map { meta, asm1, asm2 -> [meta, [asm1, asm2]] }
        .map { meta, fasta ->
            def key = groupKey(meta, fasta.size())
            [ key, fasta ]
        }
        .transpose()

    BGZIP_POLISHED(ch_bgzip_input)
    ch_versions = ch_versions.mix(BGZIP_POLISHED.out.versions)

    //
    // Logic: combine all polishing outputs into a single map for ease of publishing
    //
    ch_polishing_output = BGZIP_POLISHED.out.output
        .groupTuple(by: 0)
        .map { key, fasta -> [key.target, fasta] }
        .join(POLISHING_10X.out.bed_chunks, by: 0)
        .join(POLISHING_10X.out.longranger_bam, by: 0)
        .join(POLISHING_10X.out.longranger_bai, by: 0)
        .join(POLISHING_10X.out.longranger_summary, by: 0)
        .join(POLISHING_10X.out.merged_vcf, by: 0)
        .join(POLISHING_10X.out.merged_vcf_tbi, by: 0)
        .map { meta, fasta, bed_chunks, lr_bam, lr_bai, lr_csv, vcf, tbi ->
            return [
                id: meta.hash,
                fasta: fasta,
                bed_chunks: bed_chunks,
                longranger_bam: bam,
                longranger_bai: bai,
                longranger_csv: lr_csv,
                merged_vcf: vcf,
                merged_vcf_tbi: tbi
            ]
        }

    emit:
    polished_assemblies = ch_assemblies_polished_split
    polishing_output = ch_polishing_output
    versions = ch_versions
}













// ch_all_assemblies_after_polishing = ch_assemblies_to_polish.no_polish
//     .mix(ch_assemblies_polished)
