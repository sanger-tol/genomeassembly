def readYAML(yamlfile) {
    return new org.yaml.snakeyaml.Yaml().load(yamlfile.text)
}

def hmmify(hmmfile) {
    def extensions = [".h3f", ".h3i", ".h3m", ".h3p"]
    def hmmfiles   = extensions.collect { hmmfile + "$it" }
    return [hmmfile] + hmmfiles
}

process READ_YAML {
    executor "local"

    input:
    val(yaml)

    output:
    // Sequencing data
    tuple val(meta), val(pacbio_reads)  , emit: pacbio_reads
    tuple val(meta), val(ont_reads)     , emit: ont_reads
    tuple val(meta), val(hic_reads)     , emit: hic_reads
    tuple val(meta), val(illumina_10x)  , emit: illumina_10x_reads
    tuple val(meta), val(maternal_reads), emit: maternal_reads
    tuple val(meta), val(paternal_reads), emit: paternal_reads
    // Databases
    val(busco_lineage)                  , emit: busco_lineage
    val(oatk_mito)                      , emit: oatk_mito_hmm
    val(oatk_plastid)                   , emit: oatk_plastid_hmm

    exec:
    // Read input
    def input = readYAML(yaml)

    // Generate meta map
    meta           = input.metadata

    // Process input files
    // Sequencing data
    def sd         = input.sequencing_data
    def ill        = sd?.illumina              ? sd.illumina                      : []

    pacbio_reads   = sd?.pacbio?.reads         ? sd?.pacbio?.reads.flatten()      : []
    ont_reads      = sd?.ont?.reads            ? sd?.ont?.reads.flatten()         : []
    hic_reads      = sd?.hic?.reads            ? sd.hic.reads.flatten()           : []
    illumina_10x   = ill?.illumina_10x?.reads  ? ill.illumina_10x.reads.flatten() : []
    maternal_reads = ill?.trio?.reads_maternal ? ill.illumina_10x.reads.flatten() : []
    paternal_reads = ill?.trio?.reads_paternal ? ill.illumina_10x.reads.flatten() : []

    // Databases
    def db = input.databases

    busco_lineage = db.busco_lineage
    oatk_mito     = db?.oatk?.mito_hmm        ? hmmify(db.oatk.mito_hmm)        : []
    oatk_plastid  = db?.oatk?.plastid_hmm     ? hmmify(db.oatk.plastid_hmm)     : []

}
