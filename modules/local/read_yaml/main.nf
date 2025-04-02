def readYAML(yamlfile) {
    return new org.yaml.snakeyaml.Yaml().load(yamlfile.text)
}

def get_all_hmm(hmmfile) {
    def extensions = [".h3f", ".h3i", ".h3m", ".h3p"]
    def hmmfiles   = extensions.collect { hmmfile + "$it" }
    return [hmmfile] + hmmfiles
}

def get_fastk(map) {
    def has_fastk = map?.fastk

    if(has_fastk) {
        def fastk_files = new File(fastk.directory)
            .listFiles()
            .findAll {
                it.name =~ fastk.basename
            }

        def fk_hist = fastk_files.findAll { it =~ /\.hist$/ }
        def fk_ktab = fastk_files.findAll { it =~ /\.ktab(\.\d+)?$/ }

        if([fk_hist, fk_ktab].every { !it.isEmpty() } ){
            return [hist: fk_hist, ktab: fk_ktab]
        }
    }

    return [hist: [], ktab: []]
}

process READ_YAML {
    executor "local"

    input:
    val(yaml)

    output:
    // Raw sequencing data
    tuple val(lr_meta)  , val(long_reads), val(lr_fk.hist)  , val(lr_fk.ktab)  , emit: long_reads
    tuple val(hic_meta) , val(hic_reads) , val(hic_fk.hist) , val(hic_fk.ktab) , emit: hic_reads
    tuple val(i10x_meta), val(i10x_reads), val(i10x_fk.hist), val(i10x_fk.ktab), emit: i10x_reads
    tuple val(mat_meta) , val(mat_reads) , val(mat_fk.hist) , val(mat_fk.ktab) , emit: maternal_reads
    tuple val(pat_meta) , val(pat_reads) , val(pat_fk.hist) , val(pat_fk.ktab) , emit: paternal_reads
    // Databases
    val(busco_lineage), emit: busco_lineage
    val(oatk_mito)    , emit: oatk_mito_hmm
    val(oatk_plastid) , emit: oatk_plastid_hmm

    exec:
    // Read input
    def input = readYAML(yaml)

    // Generate meta map
    meta      = input.metadata
    lr_meta   = meta + [reads: "long", platform: input.sequencing_data.platform]
    hic_meta  = meta + [reads: "hic"]
    i10x_meta = meta + [reads: "i10x"]
    mat_meta  = meta + [reads: "mat"]
    pat_meta  = meta + [reads: "pat"]

    // Process input files
    // Sequencing data
    def sd         = input.sequencing_data

    // Raw sequencing data
    long_reads = sd.long_reads.reads
    hic_reads  = sd?.hic?.reads            ? sd.hic.reads.flatten()           : []
    i10x_reads = sd?.illumina_10x?.reads   ? ill.illumina_10x.reads.flatten() : []
    mat_reads  = sd?.trio?.maternal?.reads ? sd.trio.maternal.reads.flatten() : []
    pat_reads  = sd?.trio?.paternal?.reads ? sd.trio.paternal.reads.flatten() : []

    // Kmer databases
    def lr_fk   = get_fastk(sd.long_reads)
    def hic_fk  = get_fastk(sd.hic)
    def i10x_fk = get_fastk(sd.illumina_10x)
    def mat_fk  = get_fastk(sd.trio.maternal)
    def pat_fk  = get_fastk(sd.trio.paternal)

    // Databases
    def db = input.databases

    busco_lineage = db?.busco_lineage     ? db.busco_lineage : "auto_euk"
    oatk_mito     = db?.oatk?.mito_hmm    ? get_all_hmm(db.oatk.mito_hmm)     : []
    oatk_plastid  = db?.oatk?.plastid_hmm ? get_all_hmm(db.oatk.plastid_hmm)  : []

}
