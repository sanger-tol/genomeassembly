include { MITOHIFI_FINDMITOREFERENCE } from '../../modules/nf-core/mitohifi/findmitoreference/main'
include { MITOHIFI_MITOHIFI } from '../../modules/local/mitohifi'

workflow MITO {
    take:
    hifi_reads // channel: [ val(meta), datafile  ]
    contigs // channel: [ val(meta), datafile ] 
    mito_info // channel: [ val(meta), val(species), val(min_length), val(code), val(email) ]

    main:
    ch_versions = Channel.empty()

    mito_info.map{ meta, species, min_length, code, email -> species}.set{species}
    mito_info.map{ meta, species, min_length, code, email -> min_length}.set{min_length}
    mito_info.map{ meta, species, min_length, code, email -> code}.set{code}
    mito_info.map{ meta, species, min_length, code, email -> email}.set{email}
   // MITOHIFI_FINDMITOREFERENCE(species, email, min_length)
    MITOHIFI_FINDMITOREFERENCE("\"Deilephila porcellus\"", email, min_length)

    // TODO: concat reads
//    MITOHIFI_MITOHIFI( hifi_reads.map{meta, hifi -> [meta, hifi, []]}, 
//                       MITOHIFI_FINDMITOREFERENCE.out.fasta,
//                       MITOHIFI_FINDMITOREFERENCE.out.gb,
//                      code)    

//    MITOHIFI_MITOHIFI( contigs.map{meta, contigs -> [meta, [], contigs]}, 
//                       MITOHIFI_FINDMITOREFERENCE.out.fasta,
//                       MITOHIFI_FINDMITOREFERENCE.out.gb,
//                       code)    
    
//   MITOHIFI_MITOHIFI(  [[id:"ilDeiPorc1"], [], "/lustre/scratch124/tol/projects/darwin/users/kk16/development/mitohifi/mito-datasets/animal-ilDeiPorc1-c/ilDeiPorc1.contigs.fa"],
   MITOHIFI_MITOHIFI(  [[id:"ilDeiPorc1"], "/lustre/scratch124/tol/projects/darwin/users/kk16/development/mitohifi/mito-datasets/animal-ilDeiPorc1-r/ilDeiPorc1.reads.fa", []],
                       MITOHIFI_FINDMITOREFERENCE.out.fasta,
                       MITOHIFI_FINDMITOREFERENCE.out.gb,
                       code)    
    emit:

    fasta_ref = MITOHIFI_FINDMITOREFERENCE.out.fasta
    fasta_gb = MITOHIFI_FINDMITOREFERENCE.out.gb
    mito_ctg_fa = MITOHIFI_MITOHIFI.out.fasta
    mito_ctg_annot = MITOHIFI_MITOHIFI.out.gb ? MITOHIFI_MITOHIFI.out.gb : MITOHIFI_MITOHIFI.out.gff

//    primary_contigs = GFA_TO_FASTA_PRI.out.fasta
//    alternate_contigs = GFA_TO_FASTA_ALT.out.fasta
//    primary_hic_contigs = hifiasm_hic_on ? GFA_TO_FASTA_PRI_HIC.out.fasta : null
//    alternate_hic_contigs = hifiasm_hic_on ? GFA_TO_FASTA_ALT_HIC.out.fasta : null

    versions = ch_versions
}
