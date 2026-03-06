include { FASTA_10X_POLISHING_LONGRANGER_FREEBAYES } from '../../../subworkflows/sanger-tol/fasta_10x_polishing_longranger_freebayes'

include { CAT_CAT as CONCATENATE_ASSEMBLIES        } from '../../../modules/nf-core/cat/cat'
include { SEQKIT_GREP as SEPARATE_HAPLOTYPES       } from '../../../modules/nf-core/seqkit/grep'
include { TABIX_BGZIP as BGZIP_POLISHED            } from '../../../modules/nf-core/tabix/bgzip'

workflow POLISHING {
    take:
    ch_polishing_specs // channel: spec
    ch_assemblies      // channel: [spec, hap1, hap2]
    val_polishing_container_provided
    val_sequences_per_polishing_chunk

    main:
    ch_assemblies_polished_split = channel.empty()
    ch_polishing_output = channel.empty()

    // Exclude polishing from the CLI output if not enabled to clean
    // up onscreen output
    if(val_polishing_container_provided) {
        //
        // Logic: join all the assemblies with the purging specifications and
        // data, filter for those assemblies which are to be purged, then
        // map out the data for the purging subworkflow.
        //
        ch_polishing_input = ch_assemblies
            .combine(ch_polishing_specs)
            .filter { asm_meta, _asm1, _asm2, spec -> asm_meta.id == spec.prevID }
            .multiMap { _asm_meta, asm1, asm2, spec ->
                assemblies: [spec, [asm1, asm2]]
                polishing_reads: [spec, spec.data.polishing.reads]
            }

        //
        // Subworkflow: Run polishing with 10X reads using Longranger and Freebayes
        //
        FASTA_10X_POLISHING_LONGRANGER_FREEBAYES(
            ch_polishing_input.assemblies,
            ch_polishing_input.polishing_reads,
            val_sequences_per_polishing_chunk
        )

        //
        // Module: bgzip all polished assembly fasta
        //
        ch_bgzip_input = FASTA_10X_POLISHING_LONGRANGER_FREEBAYES.out.polished_assemblies
            .map { meta, fasta ->
                def key = groupKey(meta, fasta.size())
                [ key, fasta ]
            }
            .transpose()

        BGZIP_POLISHED(ch_bgzip_input)

        //
        // Logic: combine all polishing outputs into a single map for ease of publishing
        //
        ch_polishing_output = BGZIP_POLISHED.out.output
            .groupTuple(by: 0)
            .map { key, fasta -> [key.target, fasta] }
            .join(FASTA_10X_POLISHING_LONGRANGER_FREEBAYES.out.longranger_bam, by: 0)
            .join(FASTA_10X_POLISHING_LONGRANGER_FREEBAYES.out.longranger_bai, by: 0)
            .join(FASTA_10X_POLISHING_LONGRANGER_FREEBAYES.out.longranger_summary, by: 0)
            .join(FASTA_10X_POLISHING_LONGRANGER_FREEBAYES.out.merged_vcf, by: 0)
            .join(FASTA_10X_POLISHING_LONGRANGER_FREEBAYES.out.merged_vcf_tbi, by: 0)
            .map { spec, fasta, lr_bam, lr_bai, lr_csv, vcf, tbi ->
                return spec.subMap(["id", "stage", "data", "params", "tools"]) + [
                    output: [
                        polishing: [
                            fasta: fasta,
                            longranger_bam: lr_bam,
                            longranger_bai: lr_bai,
                            longranger_csv: lr_csv,
                            merged_vcf: vcf,
                            merged_vcf_tbi: tbi
                        ]
                    ]
                ]
            }
    }

    emit:
    polished_assemblies = ch_assemblies_polished_split
    polishing_output = ch_polishing_output
}
