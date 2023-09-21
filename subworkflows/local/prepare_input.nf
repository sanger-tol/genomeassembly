/*
Adapted from
https://github.com/mahesh-panchal/genomeassembler/blob/230fa883621e45ce3d1002eb615bd5a01a68e59b/subworkflows/local/prepare_input.nf
*/
nextflow.enable.dsl = 2

import org.yaml.snakeyaml.Yaml

include { GUNZIP as GUNZIP_PRI                       }  from '../../modules/nf-core/gunzip/main'
include { GUNZIP as GUNZIP_HAP                       }  from '../../modules/nf-core/gunzip/main'
include { SAMTOOLS_FAIDX as SAMTOOLS_FAIDX_PRIMARY   }  from '../../modules/nf-core/samtools/faidx/main'
include { SAMTOOLS_FAIDX as SAMTOOLS_FAIDX_HAPLOTIGS }  from '../../modules/nf-core/samtools/faidx/main'
include { SAMTOOLS_FAIDX as SAMTOOLS_FAIDX_MERGED    }  from '../../modules/nf-core/samtools/faidx/main'
include { FASTA_CONCAT                               }  from '../../modules/local/concat'

workflow PREPARE_INPUT {

    take:
    ch_input

    main:
    ch_versions = Channel.empty()

    Channel.of(ch_input).map { file -> readYAML( file ) }
        .set { ymlfile }
    
    ymlfile.multiMap{ data -> 
        dataset : (data.dataset ? data.dataset : []) 
        busco : (data.busco ? data.busco : [])
        hic_motif : (data.hic_motif ? data.hic_motif : [])
    }
    .set{ ch_yml_data }

    ch_yml_data.dataset.flatten()  
            .multiMap { data -> 
            id_ch : (data.id ? data.id : [])
            illumina_10X_ch : ( data.illumina_10X ? [ [id: data.id ], 
                                                       file(data.illumina_10X.reads, checkIfExists: true),
                                                       data.illumina_10X.kmer_pref ? data.illumina_10X.kmer_pref : [] ] 
                                : [] )
            pacbio_ch: ( data.pacbio ? [ [id: data.id ], 
                                          data.pacbio.reads.collect { file( it.reads, checkIfExists: true ) } ]
                        : [])
            hic_ch: ( data.HiC ? [ [id: data.id, datatype: "hic", read_group: "\'@RG\\tID:" + data.id  + "\\tPL:ILLUMINA" + "\\tSM:" + data.id + "\'" ],  
                                    data.HiC.reads.collect { file( it.reads, checkIfExists: true ) } ] 
                        : [])
        }
        .set{ dataset_ch }

    dataset_ch.hic_ch.combine(ch_yml_data.hic_motif)
                     .set{ hic_ch }

    ch_yml_data.busco.flatten()
            .map { data -> [ [ id: dataset_ch.id_ch ],
                         data.lineage_path ? file(data.lineage_path, checkIfExists: true) : [],
                         data.lineage ] }
            .set{ busco_ch }

    emit:
    hic            = hic_ch
    hifi           = dataset_ch.pacbio_ch
    illumina_10X   = dataset_ch.illumina_10X_ch
    busco          = busco_ch
    
    versions       = ch_versions.ifEmpty(null) // channel: [ versions.yml ]
}


def readYAML( yamlfile ) {
    return new Yaml().load( new FileReader( yamlfile.toString() ) )
}
