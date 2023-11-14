include { MITOHIFI_FINDMITOREFERENCE                    } from '../../modules/nf-core/mitohifi/findmitoreference/main'
include { MITOHIFI_MITOHIFI } from '../../modules/nf-core/mitohifi/mitohifi/main'

workflow ORGANELLES {
    take:
    input  // channel: [ val(meta), datafile  ]
    mito_info   // channel: [ val(species), val(min_length), val(code), val(email) ]

    main:
    ch_versions = Channel.empty()

    //
    // LOGIC: SEPARATE INPUT INTO CHANNELS
    //
    mito_info.map{ species, min_length, code, email -> species}.set{species}
    mito_info.map{ species, min_length, code, email -> min_length}.set{min_length}
    mito_info.map{ species, min_length, code, email -> code}.set{code}
    mito_info.map{ species, min_length, code, email -> email}.set{email}

    //
    // MODULE: DOWNLOAD REFERENCE ORGANELLE ASSEMBLY
    //
    MITOHIFI_FINDMITOREFERENCE(species, email, min_length)
    ch_versions = ch_versions.mix(MITOHIFI_FINDMITOREFERENCE.out.versions.first())
        
    //
    // MODULE: IDENTIFY ORGANELLE IN THE DATASET
    //
    MITOHIFI_MITOHIFI( input, 
                   MITOHIFI_FINDMITOREFERENCE.out.fasta,
                   MITOHIFI_FINDMITOREFERENCE.out.gb,
                   code)    
    ch_versions = ch_versions.mix(MITOHIFI_FINDMITOREFERENCE.out.versions.first())

    emit:

    versions = ch_versions
}
