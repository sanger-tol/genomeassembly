include { CAT_CAT as CAT_CAT_READS   } from "../../modules/nf-core/cat/cat/main"
include { FASTK_FASTK                } from "../../modules/nf-core/fastk/fastk/main"
include { FASTK_HISTEX               } from '../../modules/nf-core/fastk/histex/main'
include { GENESCOPEFK                } from "../../modules/nf-core/genescopefk/main"

workflow GENOMESCOPE_MODEL {

    take:
    reads // [meta, [reads]] 

    main: 
    ch_versions = Channel.empty()

    //
    // MODULE: MERGE ALL READS IN ONE FILE
    //
    CAT_CAT_READS( reads )
    ch_versions = ch_versions.mix(CAT_CAT_READS.out.versions)
    
    //
    // LOGIC: KEEP THE CORRECT EXTENSION
    //
    CAT_CAT_READS.out.file_out.map{ meta, reads_ch -> reads_ch.getName().endsWith('gz')
                                    ? [meta, reads_ch.getParent().toString() + '/' + reads_ch.getBaseName().toString() + '.fa.gz'] : [meta, reads_ch.getParent().toString() + '/' + reads_ch.getBaseName().toString() + '.fa'] }
                            .set{ reads_merged_ch }

    //
    // LOGIC: MAKE SURE MERGED READS HAVE THE PROPER EXTENTION
    //
    CAT_CAT_READS.out.file_out.join(reads_merged_ch)
                            .map{ meta, reads_old, reads_new -> reads_old.renameTo(reads_new); }

    //
    // MODULE: GENERATE KMER DATABASE
    //
    FASTK_FASTK( reads_merged_ch )
    ch_versions = ch_versions.mix(FASTK_FASTK.out.versions)

    //
    // MODULE: KEEP THE KMERS HISTOGRAM
    //
    FASTK_HISTEX( FASTK_FASTK.out.hist )
    ch_versions = ch_versions.mix(FASTK_HISTEX.out.versions)

    //
    // MODULE: GENERATE GENOMESCOPE KMER COVERAGE MODEL
    //
    GENESCOPEFK( FASTK_HISTEX.out.hist )
    ch_versions = ch_versions.mix(GENESCOPEFK.out.versions)

    emit:
    model = GENESCOPEFK.out.model
    hist = FASTK_FASTK.out.hist
    ktab = FASTK_FASTK.out.ktab
    
    versions = ch_versions
}

