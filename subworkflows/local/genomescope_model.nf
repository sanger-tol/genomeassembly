include { CAT_CAT as CAT_CAT_READS                  } from "../../modules/nf-core/cat/cat/main"
include { FASTK_FASTK                               } from "../../modules/nf-core/fastk/fastk/main"
include { FASTK_HISTEX                              } from '../../modules/nf-core/fastk/histex/main'
include { GENESCOPEFK                               } from "../../modules/nf-core/genescopefk/main"

workflow GENOMESCOPE_MODEL {

    take:
    reads // [meta, [reads]] 

    main: 
    reads.flatMap { meta, reads -> reads instanceof List ? reads.collect{ [ meta, it ] } : [ [ meta, reads ] ] }
                    .set{ reads_ch }

    CAT_CAT_READS( reads_ch )
    CAT_CAT_READS.out.file_out.map{ meta, path -> path.renameTo(path.getParent() + '/' + path.getBaseName() + '.fa.gz');
                                              [meta, path.getParent() + '/' + path.getBaseName() + '.fa.gz'] }
                     .set{ reads_merged_ch }
    FASTK_FASTK( reads_merged_ch )
    FASTK_HISTEX( FASTK_FASTK.out.hist )
    GENESCOPEFK ( FASTK_HISTEX.out.hist )

    emit:
    model = GENESCOPEFK.out.model
}

