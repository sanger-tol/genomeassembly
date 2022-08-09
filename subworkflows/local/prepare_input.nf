nextflow.enable.dsl = 2

import org.yaml.snakeyaml.Yaml

workflow PREPARE_INPUT {

    take:
    ch_input

    main:

    ch_input.map { file -> readYAML( file ) }
        .multiMap { file ->
            samples_ch: file.samples
        }
        .set { ymlfile }

    ymlfile.samples_ch
        .flatten()
        .dump( tag: 'YAML Samples' )
        .multiMap { data ->
            assembly_ch : ( data.assembly ? [ [ id: data.id ], file( data.assembly.primary, checkIfExists: true ), file( data.assembly.haplotigs, checkIfExists: true) ] : [] )
            illumina_ch : ( data.illumina ? [ [id: data.id ], file(data.illumina.reads, checkIfExists: true), data.illumina.kmer_pref  ] : [] )
            pacbio_ch: ( data.pacbio ? [ [id: data.id ],  data.pacbio.reads.collect { file( it.reads, checkIfExists: true ) }, data.pacbio.kmer_pref ] : [])
            hic_ch: ( data.HiC ? [ [id: data.id ],  data.HiC.reads.collect { file( it.reads, checkIfExists: true ) }, data.HiC.arima_motif ] : [])
            busco_ch : ( data.busco ? [ [id: data.id ], file(data.busco.lineages_path, checkIfExists: true), data.busco.lineage ] : [] )
        }
        .set{ yml_input }
    emit:
    assemblies = yml_input.assembly_ch.dump( tag: 'Input: Assemblies' )
    hic        = yml_input.hic_ch.dump( tag: 'Input: Hi-C' )
    hifi       = yml_input.pacbio_ch.dump( tag: 'Input: PacBio HiFi' )
    illumina   = yml_input.illumina_ch.dump( tag: 'Input: Illumina' )
    busco      = yml_input.busco_ch.dump( tag: 'Input: BUSCO lib')

}


def readYAML( yamlfile ) {
    return new Yaml().load( new FileReader( yamlfile.toString() ) )
}

/*
Notes:
Use 
```
nextflow run test.nf
```
*/
