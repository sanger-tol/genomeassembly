include { MITOHIFI_FINDMITOREFERENCE    } from '../../modules/nf-core/mitohifi/findmitoreference/main'
include { MITOHIFI_MITOHIFI as MITOHIFI_MITOHIFI_READS       } from '../../modules/nf-core/mitohifi/mitohifi/main'
include { MITOHIFI_MITOHIFI as MITOHIFI_MITOHIFI_CONTIGS     } from '../../modules/nf-core/mitohifi/mitohifi/main'

include { OATK                          } from '../../modules/nf-core/oatk/main'


workflow ORGANELLES {
    take:
    reads_input  // channel: [ val(meta), datafile  ]
    contigs_input  // channel: [ val(meta), datafile  ]
    mito_info   // channel: [ val(species), val(min_length), val(code), val(email), val(fam) ]
    plastid_info // channel: [ val(fam) ]

    main:
    ch_versions = Channel.empty()

    //
    // LOGIC: SEPARATE INPUT INTO CHANNELS
    //
    mito_info.map{ species, min_length, code, email, fam -> species}.set{species}
    mito_info.map{ species, min_length, code, email, fam -> min_length}.set{min_length}
    mito_info.map{ species, min_length, code, email, fam -> code}.set{code}
    mito_info.map{ species, min_length, code, email, fam -> email}.set{email}

    //
    // MODULE: DOWNLOAD REFERENCE ORGANELLE ASSEMBLY
    //
    MITOHIFI_FINDMITOREFERENCE(species, email, min_length)
    ch_versions = ch_versions.mix(MITOHIFI_FINDMITOREFERENCE.out.versions.first())
        
    //
    // MODULE: IDENTIFY ORGANELLE IN THE READS DATASET
    //
    MITOHIFI_MITOHIFI_READS( reads_input, 
                   MITOHIFI_FINDMITOREFERENCE.out.fasta,
                   MITOHIFI_FINDMITOREFERENCE.out.gb,
                   code)    
    ch_versions = ch_versions.mix(MITOHIFI_FINDMITOREFERENCE.out.versions.first())

    //
    // MODULE: IDENTIFY ORGANELLE IN THE ASSEMBLY DATASET
    //
    MITOHIFI_MITOHIFI_CONTIGS( contigs_input, 
                   MITOHIFI_FINDMITOREFERENCE.out.fasta,
                   MITOHIFI_FINDMITOREFERENCE.out.gb,
                   code)    
    ch_versions = ch_versions.mix(MITOHIFI_FINDMITOREFERENCE.out.versions.first())

    //
    // LOGIC: PREPARE OATK INPUT
    //
    mito_info.map{ species, min_length, code, email, fam -> [ fam ] + ['h3f', 'h3i', 'h3m', 'h3p'].collect {fam.resolveSibling(fam.name + '.' + it)} }
                                                            .set { mito_hmm_input }
    plastid_info.map{ fam -> fam ? ([ fam ] + ['h3f', 'h3i', 'h3m', 'h3p'].collect {fam.resolveSibling(fam.name + '.' + it)} ) : [[],[],[],[],[]] }
                               .set { plastid_hmm_input }
    //
    // MODULE: RUN OATK TO IDENTIFY MITO
    //
    OATK(reads_input, mito_hmm_input, plastid_hmm_input)
    ch_versions = ch_versions.mix(OATK.out.versions.first())

    emit:

    versions = ch_versions
}
