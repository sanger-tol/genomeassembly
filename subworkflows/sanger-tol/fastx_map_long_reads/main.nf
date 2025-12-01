include { FASTXALIGN_PYFASTXINDEX    } from '../../../modules/sanger-tol/fastxalign/pyfastxindex/main'
include { FASTXALIGN_MINIMAP2ALIGN   } from '../../../modules/sanger-tol/fastxalign/minimap2align/main'
include { MINIMAP2_INDEX             } from '../../../modules/nf-core/minimap2/index/main'

include { BAM_SAMTOOLS_MERGE_MARKDUP } from '../bam_samtools_merge_markdup/main'

workflow FASTX_MAP_LONG_READS {

    take:
    ch_assemblies             // Channel [meta, assembly]
    ch_fasta                  // Channel [meta, fasta] OR [meta, [fasta1, fasta2, ..., fasta_n]]
    val_reads_per_fasta_chunk // integer: Number of reads per FASTA chunk for mapping
    val_output_bam            // boolean: if true output alignments in BAM format

    main:
    ch_versions = channel.empty()

    //
    // Logic: rolling check of assembly meta objects to detect duplicates
    //
    def val_asm_meta_list = Collections.synchronizedSet(new HashSet())

    ch_assemblies
        .map { meta, _sample ->
            if (!val_asm_meta_list.add(meta)) {
                error("Error: Duplicate meta object found in `ch_assemblies` in FASTX_MAP_LONG_READS: ${meta}")
            }
            meta
        }

    //
    // Module: Index FASTA files
    //
    FASTXALIGN_PYFASTXINDEX(ch_fasta.transpose())
    ch_versions = ch_versions.mix(FASTXALIGN_PYFASTXINDEX.out.versions)

    //
    // Logic: Identify FASTA chunks
    //
    ch_fastx_chunks = FASTXALIGN_PYFASTXINDEX.out.index
        .map { meta, fasta, index, count ->
            def intcount = count.toInteger()
            def size     = val_reads_per_fasta_chunk
            def n_bins   = Math.ceil(intcount / size).toInteger()
            def chunkn   = (0..<n_bins).collect()
            def slices   = chunkn.collect { chunk ->
                def lower = chunk * size
                def upper = [lower + size, intcount].min()

                return [ lower, upper ]
            }

            return [ meta, fasta, index, chunkn, slices ]
        }

    //
    // Logic: Count the total number of cram chunks for downstream grouping
    //
    ch_n_fasta_chunks = ch_fastx_chunks
        .map { meta, _fasta, _index, chunkn, _slices -> [ meta, chunkn ] }
        .transpose()
        .groupTuple(by: 0)
        .map { meta, chunkns -> [ meta, chunkns.size() ] }

    //
    // MODULE: generate minimap2 mmi file
    //
    MINIMAP2_INDEX(ch_assemblies)
    ch_versions = ch_versions.mix(MINIMAP2_INDEX.out.versions)

    //
    // Module: Map slices of each FASTA file to the reference
    //
    ch_fasta_with_slices = ch_fastx_chunks
        .combine(MINIMAP2_INDEX.out.index, by: 0)
        .transpose()

    FASTXALIGN_MINIMAP2ALIGN(
        ch_fasta_with_slices,
        val_output_bam
    )
    ch_versions = ch_versions.mix(FASTXALIGN_MINIMAP2ALIGN.out.versions)

    //
    // Logic: Group all PAF files together, using a groupKey to output when
    //        reaching the expected count of PAF files
    //
    ch_grouped_paf = FASTXALIGN_MINIMAP2ALIGN.out.paf
        .combine(ch_n_fasta_chunks, by: 0)
        .map { meta, paf, n_chunks ->
            def key = groupKey(meta, n_chunks)
            [key, paf]
        }
        .groupTuple(by: 0, sort: { paf -> paf.getName() } )
        .map { key, paf -> [key.target, paf] } // Get meta back out of groupKey

    //
    // Logic: Group all BAM files together for merging, using a groupKey to
    //        output when reaching the expected count of PAF files
    //
    ch_merge_input = FASTXALIGN_MINIMAP2ALIGN.out.bam
        .combine(ch_n_fasta_chunks, by: 0)
        .map { meta, bam, n_chunks ->
            def key = groupKey(meta, n_chunks)
            [key, bam]
        }
        .groupTuple(by: 0, sort: { bam -> bam.getName() } )
        .map { key, bam -> [key.target, bam] } // Get meta back out of groupKey

    //
    // Logic: Wrap this in the conditional so we don't unnecessarily run
    //        samtools faidx if no bam output
    //
    ch_output_bam       = channel.empty()
    ch_output_bam_index = channel.empty()

    if(val_output_bam) {
        BAM_SAMTOOLS_MERGE_MARKDUP(
            ch_merge_input,
            ch_assemblies,
            false
        )
        ch_versions = ch_versions.mix(BAM_SAMTOOLS_MERGE_MARKDUP.out.versions)

        ch_output_bam       = ch_output_bam.mix(BAM_SAMTOOLS_MERGE_MARKDUP.out.bam)
        ch_output_bam_index = ch_output_bam.mix(BAM_SAMTOOLS_MERGE_MARKDUP.out.bam_index)
    }

    emit:
    bam       = ch_output_bam
    bam_index = ch_output_bam_index
    paf       = ch_grouped_paf
    versions  = ch_versions
}
