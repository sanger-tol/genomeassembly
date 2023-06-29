include { MITOHIFI_FINDMITOREFERENCE                    } from '../../modules/nf-core/mitohifi/findmitoreference/main'
include { MITOHIFI_MITOHIFI as MITOHIFI_MITOHIFI_READS  } from '../../modules/local/mitohifi'
include { MITOHIFI_MITOHIFI as MITOHIFI_MITOHIFI_ASM    } from '../../modules/local/mitohifi'

workflow ORGANELLES {
    take:
    hifi_reads  // channel: [ val(meta), datafile  ]
    contigs     // channel: [ val(meta), datafile ] 
    mito_info   // channel: [ val(meta), val(species), val(min_length), val(code), val(email) ]

    main:
    ch_versions = Channel.empty()

    mito_info.map{ meta, species, min_length, code, email -> species}.set{species}
    mito_info.map{ meta, species, min_length, code, email -> min_length}.set{min_length}
    mito_info.map{ meta, species, min_length, code, email -> code}.set{code}
    mito_info.map{ meta, species, min_length, code, email -> email}.set{email}
    MITOHIFI_FINDMITOREFERENCE(species, email, min_length)

    // TODO: concat reads
    if ( hifi_reads ) {
        MITOHIFI_MITOHIFI_READS( hifi_reads.map{meta, hifi -> [meta, hifi, []]}, 
                        MITOHIFI_FINDMITOREFERENCE.out.fasta,
                        MITOHIFI_FINDMITOREFERENCE.out.gb,
                        code)    
    }
    if ( contigs ) {
        MITOHIFI_MITOHIFI_ASM( contigs.map{meta, contigs -> [meta, [], contigs]}, 
                        MITOHIFI_FINDMITOREFERENCE.out.fasta,
                        MITOHIFI_FINDMITOREFERENCE.out.gb,
                        code)    
    }
    
    emit:

    fasta_ref = MITOHIFI_FINDMITOREFERENCE.out.fasta
    fasta_gb = MITOHIFI_FINDMITOREFERENCE.out.gb

    mito_ctg_fa_r = MITOHIFI_MITOHIFI_READS.out.fasta
    mito_ctg_annot_r = MITOHIFI_MITOHIFI_READS.out.gb ? MITOHIFI_MITOHIFI_READS.out.gb : MITOHIFI_MITOHIFI_READS.out.gff
    all_potential_contigs_r = MITOHIFI_MITOHIFI_READS.out.all_potential_contigs ? MITOHIFI_MITOHIFI_READS.out.all_potential_contigs: []
    contigs_annotations_r = MITOHIFI_MITOHIFI_READS.out.contigs_annotations ? MITOHIFI_MITOHIFI_READS.out.contigs_annotations : []
    contigs_circularization_r = MITOHIFI_MITOHIFI_READS.out.contigs_circularization ? MITOHIFI_MITOHIFI_READS.out.contigs_circularization : []
    contigs_filtering_r = MITOHIFI_MITOHIFI_READS.out.contigs_filtering ? MITOHIFI_MITOHIFI_READS.out.contigs_filtering : []
    coverage_mapping_r = MITOHIFI_MITOHIFI_READS.out.coverage_mapping ? MITOHIFI_MITOHIFI_READS.out.coverage_mapping : []
    coverage_plot_r = MITOHIFI_MITOHIFI_READS.out.coverage_plot ? MITOHIFI_MITOHIFI_READS.out.coverage_plot : []
    final_mitogenome_annotation_r = MITOHIFI_MITOHIFI_READS.out.final_mitogenome_annotation ? MITOHIFI_MITOHIFI_READS.out.final_mitogenome_annotation : []
    final_mitogenome_choice_r = MITOHIFI_MITOHIFI_READS.out.final_mitogenome_choice ? MITOHIFI_MITOHIFI_READS.out.final_mitogenome_choice : []
    final_mitogenom_coverage_r = MITOHIFI_MITOHIFI_READS.out.final_mitogenome_coverage ? MITOHIFI_MITOHIFI_READS.out.final_mitogenome_coverage : []
    potential_contigs_r = MITOHIFI_MITOHIFI_READS.out.potential_contigs ? MITOHIFI_MITOHIFI_READS.out.potential_contigs : []
    reads_mapping_and_assembly_r = MITOHIFI_MITOHIFI_READS.out.reads_mapping_and_assembly ? MITOHIFI_MITOHIFI_READS.out.reads_mapping_and_assembly : []
    shared_genes_r = MITOHIFI_MITOHIFI_READS.out.shared_genes ? MITOHIFI_MITOHIFI_READS.out.shared_genes : []
    
    mito_ctg_fa_asm = MITOHIFI_MITOHIFI_ASM.out.fasta
    mito_ctg_annot_asm = MITOHIFI_MITOHIFI_ASM.out.gb ? MITOHIFI_MITOHIFI_ASM.out.gb : MITOHIFI_MITOHIFI_ASM.out.gff
    all_potential_contigs_asm = MITOHIFI_MITOHIFI_ASM.out.all_potential_contigs ? MITOHIFI_MITOHIFI_ASM.out.all_potential_contigs: []
    contigs_annotations_asm = MITOHIFI_MITOHIFI_ASM.out.contigs_annotations ? MITOHIFI_MITOHIFI_ASM.out.contigs_annotations : []
    contigs_circularization_asm = MITOHIFI_MITOHIFI_ASM.out.contigs_circularization ? MITOHIFI_MITOHIFI_ASM.out.contigs_circularization : []
    contigs_filtering_asm = MITOHIFI_MITOHIFI_ASM.out.contigs_filtering ? MITOHIFI_MITOHIFI_ASM.out.contigs_filtering : []
    coverage_mapping_asm = MITOHIFI_MITOHIFI_ASM.out.coverage_mapping ? MITOHIFI_MITOHIFI_ASM.out.coverage_mapping : []
    coverage_plot_asm = MITOHIFI_MITOHIFI_ASM.out.coverage_plot ? MITOHIFI_MITOHIFI_ASM.out.coverage_plot : []
    final_mitogenome_annotation_asm = MITOHIFI_MITOHIFI_ASM.out.final_mitogenome_annotation ? MITOHIFI_MITOHIFI_ASM.out.final_mitogenome_annotation : []
    final_mitogenome_choice_asm = MITOHIFI_MITOHIFI_ASM.out.final_mitogenome_choice ? MITOHIFI_MITOHIFI_ASM.out.final_mitogenome_choice : []
    final_mitogenom_coverage_asm = MITOHIFI_MITOHIFI_ASM.out.final_mitogenome_coverage ? MITOHIFI_MITOHIFI_ASM.out.final_mitogenome_coverage : []
    potential_contigs_asm = MITOHIFI_MITOHIFI_ASM.out.potential_contigs ? MITOHIFI_MITOHIFI_ASM.out.potential_contigs : []
    reads_mapping_and_assembly_asm = MITOHIFI_MITOHIFI_ASM.out.reads_mapping_and_assembly ? MITOHIFI_MITOHIFI_ASM.out.reads_mapping_and_assembly : []
    shared_genes_asm = MITOHIFI_MITOHIFI_ASM.out.shared_genes ? MITOHIFI_MITOHIFI_ASM.out.shared_genes : []

    versions = ch_versions
}
