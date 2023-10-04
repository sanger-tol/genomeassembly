include { HIFIASM as HIFIASM_PRI                } from '../../modules/nf-core/hifiasm/main'
include { HIFIASM as HIFIASM_HIC                } from '../../modules/nf-core/hifiasm/main'

include { GFA_TO_FASTA as GFA_TO_FASTA_PRI      } from '../../modules/local/gfa_to_fasta'
include { GFA_TO_FASTA as GFA_TO_FASTA_ALT      } from '../../modules/local/gfa_to_fasta'
include { GFA_TO_FASTA as GFA_TO_FASTA_PRI_HIC  } from '../../modules/local/gfa_to_fasta'
include { GFA_TO_FASTA as GFA_TO_FASTA_ALT_HIC  } from '../../modules/local/gfa_to_fasta'

workflow RAW_ASSEMBLY {
    take:
    hifi_reads // channel: [ val(meta), [ datafile ] ]
    hic_reads // channel: [ datafile ] 
    hifiasm_hic_on // val: True/False

    main:
    ch_versions = Channel.empty()

    HIFIASM_PRI(hifi_reads, [], [], [], [], [])
    ch_versions = ch_versions.mix(HIFIASM_PRI.out.versions)

    GFA_TO_FASTA_PRI( HIFIASM_PRI.out.primary_contigs )
    GFA_TO_FASTA_ALT( HIFIASM_PRI.out.alternate_contigs )
    ch_versions = ch_versions.mix(GFA_TO_FASTA_PRI.out.versions)

    if ( hifiasm_hic_on ) {
        HIFIASM_HIC(hifi_reads, [], [], [], [], hic_reads)
        GFA_TO_FASTA_PRI_HIC( HIFIASM_HIC.out.hic_primary_contigs )
        GFA_TO_FASTA_ALT_HIC( HIFIASM_HIC.out.hic_alternate_contigs )
    }

    emit:
    raw_unitigs = HIFIASM_PRI.out.raw_unitigs
    source_overlaps = HIFIASM_PRI.out.source_overlaps   
    reverse_overlaps = HIFIASM_PRI.out.reverse_overlaps
    corrected_reads = HIFIASM_PRI.out.corrected_reads
    primary_contigs_gfa = HIFIASM_PRI.out.primary_contigs
    alternate_contigs_gfa = HIFIASM_PRI.out.alternate_contigs
    processed_unitigs = HIFIASM_PRI.out.processed_unitigs
    
    primary_hic_contigs_gfa = hifiasm_hic_on ? HIFIASM_HIC.out.hic_primary_contigs : null
    alternate_hic_contigs_gfa = hifiasm_hic_on ? HIFIASM_HIC.out.hic_alternate_contigs : null 
    phased_hic_contigs_hap1_gfa = hifiasm_hic_on ? HIFIASM_HIC.out.paternal_contigs : null
    phased_hic_contigs_hap2_gfa = hifiasm_hic_on ? HIFIASM_HIC.out.maternal_contigs : null

    primary_contigs = GFA_TO_FASTA_PRI.out.fasta
    alternate_contigs = GFA_TO_FASTA_ALT.out.fasta
    primary_hic_contigs = hifiasm_hic_on ? GFA_TO_FASTA_PRI_HIC.out.fasta : null
    alternate_hic_contigs = hifiasm_hic_on ? GFA_TO_FASTA_ALT_HIC.out.fasta : null

    versions = ch_versions
}
