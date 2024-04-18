/*
Adapted from
https://github.com/mahesh-panchal/genomeassembler/blob/230fa883621e45ce3d1002eb615bd5a01a68e59b/subworkflows/local/prepare_input.nf
*/
nextflow.enable.dsl = 2

import org.yaml.snakeyaml.Yaml

workflow PREPARE_INPUT {

    take:
    ch_input

    main:
    ch_versions = Channel.empty()

    //
    // LOAD YAML
    //
    Channel.of(ch_input).map { file -> readYAML( file ) }
        .set { ymlfile }
    
    //
    // LOGIC: DIVIDE INPUT INTO BLOCKS BY SEMANTICS
    //
    ymlfile.multiMap{ data -> 
        dataset : (data.dataset ? data.dataset : []) 
        busco : (data.busco ? data.busco : [])
        mito: ( data.mito ? ['\"'+data.mito.species+'\"', data.mito.min_length, data.mito.code, data.mito.email ? data.mito.email : "\"\"", data.mito.fam ? data.mito.fam : "\"\"" ] : [])
        plastid : ( data.plastid ? ( data.plastid.fam ? data.plastid.fam : "\"\"" ) : [])
        hic_motif : (data.hic_motif ? data.hic_motif : [])
    }
    .set{ ch_yml_data }

    //
    // LOGIC: DIVIDE DATASET INTO BLOCKS BY DATATYPE
    //
    ch_yml_data.dataset.flatten()  
            .multiMap { data -> 
            id_ch : (data.id ? [id: data.id] : [])
            illumina_10X_ch : ( data.illumina_10X ? [ [id: data.id ], 
                                                       file(data.illumina_10X.reads, checkIfExists: true),
                                                       data.illumina_10X.kmer_pref ? data.illumina_10X.kmer_pref : [] ] 
                                : [] )
            pacbio_ch: ( data.pacbio ? [ [id: data.id ], 
                                          data.pacbio.reads.collect { file( it.reads, checkIfExists: true ) } ]
                        : [])
            hic_ch: ( data.HiC ? [ [id: data.id ],  
                                    data.HiC.reads.collect { file( it.reads, checkIfExists: true ) } ] 
                        : [])
        }
        .set{ dataset_ch }

    //
    // LOGIC: ADD HIC MOTIF TO DATASET HIC CHANNEL
    //
    dataset_ch.hic_ch.combine(ch_yml_data.hic_motif)
                     .set{ hic_ch }

    //
    // LOGIC: REFACTOR BUSCO CHANNEL TO ADD META
    //
    dataset_ch.id_ch.combine (
        ch_yml_data.busco.flatten()
                .map { data -> [
                             data.lineage_path ? file(data.lineage_path, checkIfExists: true) : [],
                             data.lineage ] }
        )
            .set{ busco_ch }

    emit:

    hic            = hic_ch
    hifi           = dataset_ch.pacbio_ch
    illumina_10X   = dataset_ch.illumina_10X_ch
    busco          = busco_ch
    mito           = ch_yml_data.mito 
    plastid        = ch_yml_data.plastid
    versions       = ch_versions.ifEmpty(null) // channel: [ versions.yml ]
}


def readYAML( yamlfile ) {
    return new Yaml().load( new FileReader( yamlfile.toString() ) )
}
