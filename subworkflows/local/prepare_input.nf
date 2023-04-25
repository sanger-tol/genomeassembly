/*
Adapted from
https://github.com/mahesh-panchal/genomeassembler/blob/230fa883621e45ce3d1002eb615bd5a01a68e59b/subworkflows/local/prepare_input.nf
*/
nextflow.enable.dsl = 2

import org.yaml.snakeyaml.Yaml

include { GUNZIP as GUNZIP_PRI }  from '../../modules/nf-core/gunzip/main'
include { GUNZIP as GUNZIP_HAP }  from '../../modules/nf-core/gunzip/main'
include { SAMTOOLS_FAIDX as SAMTOOLS_FAIDX_PRIMARY   }  from '../../modules/nf-core/samtools/faidx/main'
include { SAMTOOLS_FAIDX as SAMTOOLS_FAIDX_HAPLOTIGS }  from '../../modules/nf-core/samtools/faidx/main'
include { SAMTOOLS_FAIDX as SAMTOOLS_FAIDX_MERGED    }  from '../../modules/nf-core/samtools/faidx/main'
include { FASTA_CONCAT }  from '../../modules/local/concat'

workflow PREPARE_INPUT {

    take:
    ch_input

    main:
    ch_versions = Channel.empty()

    ch_input.map { file -> readYAML( file ) }
        .set { ymlfile }

    ymlfile.flatten()
        .multiMap { data -> 
                primary_ch: ( data.assembly ? data.assembly.primary : [] )
                haplotigs_ch: ( data.assembly ? data.assembly.haplotigs ? data.assembly.haplotigs : [] : [] )
                id_ch : (data.id ? [[:] ,data.id] : [])
            }
        .set{ assembly_input }

    // Prepare primary  
    assembly_input.primary_ch.map { fasta ->
    effect = fasta.endsWith('.gz') ? 'gunzip' : 'none'
    [ ['effect':effect], fasta ]
    }
    .branch {
        meta, fasta ->
        gzip : meta.effect == "gunzip"
            return [ [:], fasta ]
        geno : meta.effect == "none"
            return [ [:], fasta ]
    }
    .set { ch_asm }

    GUNZIP_PRI ( ch_asm.gzip ).gunzip
    .concat ( ch_asm.geno )
    .flatten()
    .toList()
    .join( assembly_input.id_ch )
    .map{ blank, p, meta -> [ [id: meta], p]}
    .set { ch_asm_pri_fasta }
    ch_versions = ch_versions.mix(GUNZIP_PRI.out.versions)

    // Prepare haplotigs 
    assembly_input.haplotigs_ch.map { fasta ->
    effect = fasta.endsWith('.gz') ? 'gunzip' : 'none'
    [ ['effect':effect], fasta ]
    }
    .branch {
        meta, fasta ->
        gzip : meta.effect == "gunzip"
            return [ [:], fasta ]
        geno : meta.effect == "none"
            return [ [:], fasta ]
    }
    .set { ch_asm_hap }

    GUNZIP_HAP ( ch_asm_hap.gzip ).gunzip
    .concat ( ch_asm_hap.geno )
    .flatten()
    .toList()
    .join( assembly_input.id_ch )
    .map{ blank, p, meta -> [ [id: meta], p]}
    .set { ch_asm_hap_fasta }
    ch_versions = ch_versions.mix(GUNZIP_HAP.out.versions)

    // Index all fasta files
    SAMTOOLS_FAIDX_PRIMARY( ch_asm_pri_fasta )
    ch_versions = ch_versions.mix(SAMTOOLS_FAIDX_PRIMARY.out.versions)    
    SAMTOOLS_FAIDX_HAPLOTIGS( ch_asm_hap_fasta )
    ch_versions = ch_versions.mix(SAMTOOLS_FAIDX_HAPLOTIGS.out.versions)    


    ch_asm_pri_fasta.join(SAMTOOLS_FAIDX_PRIMARY.out.fai)
                              .set{ primary_asm }
    ch_asm_hap_fasta.join(SAMTOOLS_FAIDX_HAPLOTIGS.out.fai)
                              .set{ haplotigs_asm } 

    ymlfile.flatten()
            .multiMap { data -> 
            illumina_10X_ch : ( data.illumina_10X ? [ [id: data.id ], 
                                                       file(data.illumina_10X.reads, checkIfExists: true),
                                                       data.illumina_10X.kmer_pref ? data.illumina_10X.kmer_pref : [] ] 
                                : [] )
            pacbio_ch: ( data.pacbio ? [ [id: data.id ], 
                                          data.pacbio.reads.collect { file( it.reads, checkIfExists: true ) } ]
                        : [])
            hic_ch: ( data.HiC ? [ [id: data.id, datatype: "hic", read_group: "\'@RG\\tID:" + data.id  + "\\tPL:ILLUMINA" + "\\tSM:" + data.id + "\'" ],  
                                    data.HiC.reads.collect { file( it.reads, checkIfExists: true ) }, 
                                    data.HiC.arima_motif ] 
                        : [])
            busco_ch : ( data.busco ? [ [id: data.id ], 
                                         file(data.busco.lineages_path, checkIfExists: true),
                                         data.busco.lineage ] 
                        : [] )
        }
        .set{ yml_input }

    emit:
    primary_asm    = primary_asm
    haplotigs_asm  = haplotigs_asm
    hic            = yml_input.hic_ch
    hifi           = yml_input.pacbio_ch
    illumina_10X   = yml_input.illumina_10X_ch
    busco          = yml_input.busco_ch
    
    versions       = ch_versions.ifEmpty(null) // channel: [ versions.yml ]

}


def readYAML( yamlfile ) {
    return new Yaml().load( new FileReader( yamlfile.toString() ) ).samples
}
