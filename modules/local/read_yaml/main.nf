def readYAML(yamlfile) {
    return new org.yaml.snakeyaml.Yaml().load(yamlfile.text)
}

def get_all_hmm(hmmfile) {
    def extensions = [".h3f", ".h3i", ".h3m", ".h3p"]
    def hmmfiles   = extensions.collect { hmmfile + "$it" }
    return [hmmfile] + hmmfiles
}

def get_seq_data(map) {
    def reads    = map.reads.flatten()
    def coverage = map?.coverage ? map.coverage : -1

    def output = [
        reads: reads,
        coverage: coverage,
        hist: [],
        ktab: [],
        kmer_size: -1
    ]

    def fastk = map?.fastk
    if(fastk) {
        def fastk_files = new File(fastk.directory)
            .listFiles()
            .findAll {
                it.name =~ fastk.basename
            }

        def fk_hist = fastk_files.findAll { it =~ /\.hist$/ }
        def fk_ktab = fastk_files.findAll { it =~ /\.ktab(\.\d+)?$/ }

        if([fk_hist, fk_ktab].every { !it.isEmpty() } ){
            output = output + [hist: fk_hist, ktab: fk_ktab, kmer_size: fastk.kmer_size]
        }
    }

    return output
}

process READ_YAML {
    executor "local"

    input:
    val(yaml)

    output:
    // Raw sequencing data
    tuple val(lr_meta)  , val(lr.reads)  , val(lr.hist)  , val(lr.ktab)  , emit: long_reads
    tuple val(hic_meta) , val(hic.reads) , val(hic.hist) , val(hic.ktab) , emit: hic_reads
    tuple val(i10x_meta), val(i10x.reads), val(i10x.hist), val(i10x.ktab), emit: i10x_reads
    tuple val(mat_meta) , val(mat.reads) , val(mat.hist) , val(mat.ktab) , emit: maternal_reads
    tuple val(pat_meta) , val(pat.reads) , val(pat.hist) , val(pat.ktab) , emit: paternal_reads
    // Databases
    val(busco_lineage), emit: busco_lineage
    val(oatk_mito)    , emit: oatk_mito_hmm
    val(oatk_plastid) , emit: oatk_plastid_hmm

    exec:
    // Read input
    def input = readYAML(yaml)
    def sd    = input.sequencing_data
    def db    = input.databases

    // Process input files
    def lr   = get_seq_data(sd.long_reads)
    def hic  = sd.hic           ? get_seq_data(sd.hic)           : [reads: [], hist: [], ktab: [], kmer_size: null]
    def i10x = sd.illumina_10x  ? get_seq_data(sd.illumina_10x)  : [reads: [], hist: [], ktab: [], kmer_size: null]
    def mat  = sd.trio.maternal ? get_seq_data(sd.trio.maternal) : [reads: [], hist: [], ktab: [], kmer_size: null]
    def pat  = sd.trio.paternal ? get_seq_data(sd.trio.paternal) : [reads: [], hist: [], ktab: [], kmer_size: null]

    // Generate meta maps
    meta      = input.metadata
    lr_meta   = meta + [
        read_type: "long",
        kmer_size: lr.kmer_size,
        coverage : lr.coverage,
        platform : sd.long_reads.platform
    ]
    hic_meta  = meta + [read_type: "hic" , kmer_size: hic.kmer_size ]
    i10x_meta = meta + [read_type: "i10x", kmer_size: i10x.kmer_size]
    mat_meta  = meta + [read_type: "mat" , kmer_size: mat.kmer_size ]
    pat_meta  = meta + [read_type: "pat" , kmer_size: pat.kmer_size ]

    // Databases
    busco_lineage = db?.busco_lineage     ? db.busco_lineage : "auto_euk"
    oatk_mito     = db?.oatk?.mito_hmm    ? get_all_hmm(db.oatk.mito_hmm)     : []
    oatk_plastid  = db?.oatk?.plastid_hmm ? get_all_hmm(db.oatk.plastid_hmm)  : []

}
