/*
Adapted from
https://github.com/mahesh-panchal/genomeassembler/blob/230fa883621e45ce3d1002eb615bd5a01a68e59b/subworkflows/local/prepare_input.nf
*/
nextflow.enable.dsl = 2

import org.yaml.snakeyaml.Yaml

include { GUNZIP }  from '../../modules/nf-core/modules/gunzip/main'
include { SAMTOOLS_FAIDX as SAMTOOLS_FAIDX_PRIMARY   }  from '../../modules/nf-core/modules/samtools/faidx/main'
include { SAMTOOLS_FAIDX as SAMTOOLS_FAIDX_HAPLOTIGS }  from '../../modules/nf-core/modules/samtools/faidx/main'
include { SAMTOOLS_FAIDX as SAMTOOLS_FAIDX_MERGED    }  from '../../modules/nf-core/modules/samtools/faidx/main'
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
  
    assembly_input.primary_ch.concat( assembly_input.haplotigs_ch )
    .map { fasta ->
    todo = fasta.endsWith('.gz') ? 'gunzip' : 'none'
    [ ['todo':todo], fasta ]
    }
    .branch {
        meta, fasta ->
        gzip : meta.todo == "gunzip"
            return [ [:], fasta ]
        geno : meta.todo == "none"
            return [ [:], fasta ]
    }
    .set { ch_asm }

    GUNZIP ( ch_asm.gzip ).gunzip
    .concat ( ch_asm.geno )
    .flatten()
    .toList()
    .map{ blank1, p, blank2, h -> [ blank1, p, h ] } 
    .join( assembly_input.id_ch )
    .map{ blank, p, h, meta -> [ [id: meta], p, h]}
    .set { ch_fasta }
    ch_versions = ch_versions.mix(GUNZIP.out.versions)

    //  Merge primary and haplotigs for polishing
    FASTA_CONCAT ( ch_fasta )
    ch_fasta.join(FASTA_CONCAT.out.fasta_merged)
            .set{ ch_fasta_merged }
    ch_versions = ch_versions.mix(FASTA_CONCAT.out.versions)

    ch_fasta_merged.multiMap{ it -> 
            primary: [it[0], it[1]]
            haplotigs: [it[0], it[2]]
            merged: [it[0], it[3]]}
            .set{ch_fasta_separated}

    // Index all fasta files
    SAMTOOLS_FAIDX_PRIMARY{ ch_fasta_separated.primary }
    ch_versions = ch_versions.mix(SAMTOOLS_FAIDX_PRIMARY.out.versions)    
    SAMTOOLS_FAIDX_HAPLOTIGS(ch_fasta_separated.haplotigs)
    ch_versions = ch_versions.mix(SAMTOOLS_FAIDX_HAPLOTIGS.out.versions)    
    SAMTOOLS_FAIDX_MERGED(ch_fasta_separated.merged)
    ch_versions = ch_versions.mix(SAMTOOLS_FAIDX_MERGED.out.versions)    

    SAMTOOLS_FAIDX_PRIMARY.out.fai.join(SAMTOOLS_FAIDX_HAPLOTIGS.out.fai)
                                  .join(SAMTOOLS_FAIDX_MERGED.out.fai)
                                  .set{ ch_fai }
    ymlfile.flatten()
            .multiMap { data -> 
            illumina_10X_ch : ( data.illumina_10X ? [ [id: data.id ], 
                                                       file(data.illumina_10X.reads, checkIfExists: true),
                                                       data.illumina_10X.kmer_pref ? data.illumina_10X.kmer_pref : [] ] 
                                : [] )
            pacbio_ch: ( data.pacbio ? [ [id: data.id ], 
                                          data.pacbio.reads.collect { file( it.reads, checkIfExists: true ) }, 
                                          data.pacbio.kmer_pref ? data.pacbio.kmer_pref : [] ] 
                        : [])
            hic_ch: ( data.HiC ? [ [id: data.id, datatype: "hic", read_group: "\'@RG\\tID:" + data.id + "\\tPL:ILLUMINA" + "\\tSM:" + data.id+"\'" ],  
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
    assemblies     = ch_fasta_merged
    indices        = ch_fai
    hic            = yml_input.hic_ch
    hifi           = yml_input.pacbio_ch
    illumina_10X   = yml_input.illumina_10X_ch
    busco          = yml_input.busco_ch
    
    versions       = ch_versions.ifEmpty(null) // channel: [ versions.yml ]

}


def readYAML( yamlfile ) {
    return new Yaml().load( new FileReader( yamlfile.toString() ) ).samples
}
