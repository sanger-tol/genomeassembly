def readYAML(yamlfile) {
    return new org.yaml.snakeyaml.Yaml().load(yamlfile.text)
}

def get_all_hmm(hmmfile) {
    def extensions = [".h3f", ".h3i", ".h3m", ".h3p"]
    def hmmfiles   = extensions.collect { hmmfile + "$it" }
    return [hmmfile] + hmmfiles
}

def get_fastk_files(fastk_dir, fastk_basename) {
    def fastk_files = new File(fastk_dir)
        .listFiles()
        .findAll {
            it.name =~ fastk_basename
        }

    def hist = fastk_files.findAll { it =~ /\.hist$/ }
    def ktab = fastk_files.findAll { it =~ /\.ktab(\.\d+)?$/ }

    return [hist, ktab]
}

process READ_YAML {
    executor "local"

    input:
    val(yaml)

    output:
    // Raw sequencing data
    tuple val(pacbio_meta)  , val(pacbio_reads)                       , emit: pacbio_reads
    tuple val(ont_meta)     , val(ont_reads)                          , emit: ont_reads
    tuple val(illumina_meta), val(hic_reads)                          , emit: hic_reads
    tuple val(illumina_meta), val(illumina_10x)                       , emit: illumina_10x_reads
    tuple val(illumina_meta), val(maternal_reads)                     , emit: maternal_reads
    tuple val(illumina_meta), val(paternal_reads)                     , emit: paternal_reads
    // Kmer databases
    tuple val(pacbio_meta)  , val(pacbio_fk_hist), val(pacbio_fk_ktab), emit: pacbio_fk
    tuple val(ont_meta)     , val(ont_fk_hist)   , val(ont_fk_ktab)   , emit: ont_fk
    // Databases
    val(busco_lineage)                                                , emit: busco_lineage
    val(oatk_mito)                                                    , emit: oatk_mito_hmm
    val(oatk_plastid)                                                 , emit: oatk_plastid_hmm

    exec:
    // Read input
    def input = readYAML(yaml)

    // Generate meta map
    meta           = input.metadata
    pacbio_meta    = meta + [ platform: "pacbio_hifi" ]
    ont_meta       = meta + [ platform: "ont"         ]
    illumina_meta  = meta + [ platform: "illumina"    ]

    // Process input files
    // Sequencing data
    def sd         = input.sequencing_data
    def ill        = sd?.illumina              ? sd.illumina                        : []

    // Raw sequencing data
    pacbio_reads   = sd?.pacbio?.reads         ? sd?.pacbio?.reads.flatten()        : []
    ont_reads      = sd?.ont?.reads            ? sd?.ont?.reads.flatten()           : []
    hic_reads      = sd?.hic?.reads            ? sd.hic.reads.flatten()             : []
    illumina_10x   = ill?.illumina_10x?.reads  ? ill.illumina_10x.reads.flatten()   : []
    maternal_reads = ill?.trio?.reads_maternal ? ill.trio?.reads_maternal.flatten() : []
    paternal_reads = ill?.trio?.reads_paternal ? ill.trio?.reads_paternal.flatten() : []

    // Kmer databases
    def pacbio_fastk = sd?.pacbio?.fastk?.directory_path ?
        get_fastk_files(sd.pacbio.fastk.directory_path, sd.pacbio.fastk.basename)
        : []
    def ont_fastk  = sd?.ont?.fastk?.directory_path ?
        get_fastk_files(sd.ont.fastk.directory_path, sd.ont.fastk.basename)
        : []

    pacbio_fk_hist = pacbio_fastk              ? pacbio_fastk[0]                    : []
    pacbio_fk_ktab = pacbio_fastk              ? pacbio_fastk[1]                    : []
    ont_fk_hist    = ont_fastk                 ? ont_fastk[0]                       : []
    ont_fk_ktab    = ont_fastk                 ? ont_fastk[1]                       : []

    // Databases
    def db = input.databases

    busco_lineage = db.busco_lineage
    oatk_mito     = db?.oatk?.mito_hmm          ? get_all_hmm(db.oatk.mito_hmm)     : []
    oatk_plastid  = db?.oatk?.plastid_hmm       ? get_all_hmm(db.oatk.plastid_hmm)  : []

}
