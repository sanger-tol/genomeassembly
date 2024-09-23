include { HIFIASM as HIFIASM_PRI                } from '../../modules/nf-core/hifiasm/main'
include { HIFIASM as HIFIASM_HIC                } from '../../modules/nf-core/hifiasm/main'

include { GFA_TO_FASTA as GFA_TO_FASTA_PRI      } from '../../modules/local/gfa_to_fasta'
include { GFA_TO_FASTA as GFA_TO_FASTA_ALT      } from '../../modules/local/gfa_to_fasta'
include { GFA_TO_FASTA as GFA_TO_FASTA_HAP1_HIC  } from '../../modules/local/gfa_to_fasta'
include { GFA_TO_FASTA as GFA_TO_FASTA_HAP2_HIC  } from '../../modules/local/gfa_to_fasta'

workflow RAW_ASSEMBLY {
    take:
    hifi_reads // channel: [ val(meta), [ datafile ] ]
    hic_reads // channel: [ datafile ] 
    hifiasm_hic_on // val: True/False
    matdb
    patdb

    main:
    ch_versions = Channel.empty()

    // 
    // MODULE: RUN HIFIASM IN STANDARD WAY
    //
    HIFIASM_PRI(hifi_reads, patdb, matdb, [], [], [])
    ch_versions = ch_versions.mix(HIFIASM_PRI.out.versions)

    //
    // MODULE: CONVERT PRIMARY CONTIGS TO FASTA
    //
    GFA_TO_FASTA_PRI( HIFIASM_PRI.out.primary_contigs )

    //
    // MODULE: CONVERT ALT CONTIGS TO FASTA
    //
    GFA_TO_FASTA_ALT( HIFIASM_PRI.out.alternate_contigs )
    ch_versions = ch_versions.mix(GFA_TO_FASTA_PRI.out.versions)

    //
    // LOGIC: IF FLAG SWITCHED ON RUN HIFIASM IN HIC MODE
    //
    if ( hifiasm_hic_on ) {
        //
        // MODULE: RUN HIFIASM IN HIC MODE
        //
        HIFIASM_HIC(hifi_reads, patdb, matdb, [], [], hic_reads)
        
        //
        // MODULE: CONVERT HIFIASM-HIC PRIMARY CONTIGS TO FASTA
        //
        GFA_TO_FASTA_HAP1_HIC( HIFIASM_HIC.out.hap1_contigs )
        
        //
        // MODULE: CONVERT HIFIASM-HIC ALT CONTIGS TO FASTA
        //
        GFA_TO_FASTA_HAP2_HIC( HIFIASM_HIC.out.hap2_contigs )
    }

    emit:
    primary_contigs = GFA_TO_FASTA_PRI.out.fasta
    alternate_contigs = GFA_TO_FASTA_ALT.out.fasta
    hap1_hic_contigs = hifiasm_hic_on ? GFA_TO_FASTA_HAP1_HIC.out.fasta : null
    hap2_hic_contigs = hifiasm_hic_on ? GFA_TO_FASTA_HAP2_HIC.out.fasta : null

    versions = ch_versions
}
