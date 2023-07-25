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
        mitohifi_reads = MITOHIFI_MITOHIFI_READS.out
    }
    if ( contigs ) {
        MITOHIFI_MITOHIFI_ASM( contigs.map{meta, contigs -> [meta, [], contigs]}, 
                        MITOHIFI_FINDMITOREFERENCE.out.fasta,
                        MITOHIFI_FINDMITOREFERENCE.out.gb,
                        code)    
        mitohifi_asm = MITOHIFI_MITOHIFI_ASM.out

    }
    
    emit:

    fasta_ref = MITOHIFI_FINDMITOREFERENCE.out.fasta
    fasta_gb = MITOHIFI_FINDMITOREFERENCE.out.gb

    mito_ctg_fa_r =  mitohifi_reads ? mitohifi_reads.fasta : []
    mito_ctg_annot_r = mitohifi_reads ? mitohifi_reads.gb ? mitohifi_reads.gb : mitohifi_reads.gff ? mitohifi_reads.gff : [] : []
    all_potential_contigs_r = mitohifi_reads ? mitohifi_reads.all_potential_contigs: []
    contigs_annotations_r = mitohifi_reads ? mitohifi_reads.contigs_annotations : []
    contigs_circularization_r = mitohifi_reads ? mitohifi_reads.contigs_circularization : []
    contigs_filtering_r = mitohifi_reads ? mitohifi_reads.contigs_filtering : []
    coverage_mapping_r = mitohifi_reads ? mitohifi_reads.coverage_mapping : []
    coverage_plot_r = mitohifi_reads ? mitohifi_reads.coverage_plot : []
    final_mitogenome_annotation_r = mitohifi_reads ? mitohifi_reads.final_mitogenome_annotation : []
    final_mitogenome_choice_r = mitohifi_reads ? mitohifi_reads.final_mitogenome_choice : []
    final_mitogenom_coverage_r = mitohifi_reads ? mitohifi_reads.final_mitogenome_coverage : []
    potential_contigs_r = mitohifi_reads ? mitohifi_reads.potential_contigs : []
    reads_mapping_and_assembly_r = mitohifi_reads ? mitohifi_reads.reads_mapping_and_assembly : []
    shared_genes_r = mitohifi_reads ? mitohifi_reads.shared_genes : []
    
    mito_ctg_fa_asm = mitohifi_asm ? mitohifi_asm.fasta : []
    mito_ctg_annot_asm = mitohifi_asm ? mitohifi_asm.gb ? mitohifi_asm.gb : mitohifi_asm.gff ? mitohifi_asm.gff : [] : []
    all_potential_contigs_asm = mitohifi_asm ? mitohifi_asm.all_potential_contigs: []
    contigs_annotations_asm = mitohifi_asm ? mitohifi_asm.contigs_annotations : []
    contigs_circularization_asm = mitohifi_asm ? mitohifi_asm.contigs_circularization : []
    contigs_filtering_asm = mitohifi_asm ? mitohifi_asm.contigs_filtering : []
    coverage_mapping_asm = mitohifi_asm ? mitohifi_asm.coverage_mapping : []
    coverage_plot_asm = mitohifi_asm ? mitohifi_asm.coverage_plot : []
    final_mitogenome_annotation_asm = mitohifi_asm ? mitohifi_asm.final_mitogenome_annotation : []
    final_mitogenome_choice_asm = mitohifi_asm ? mitohifi_asm.final_mitogenome_choice : []
    final_mitogenom_coverage_asm = mitohifi_asm ? mitohifi_asm.final_mitogenome_coverage : []
    potential_contigs_asm = mitohifi_asm ? mitohifi_asm.potential_contigs : []
    reads_mapping_and_assembly_asm = mitohifi_asm ? mitohifi_asm.reads_mapping_and_assembly : []
    shared_genes_asm = mitohifi_asm.shared_genes ? mitohifi_asm.shared_genes : []

    versions = ch_versions
}
