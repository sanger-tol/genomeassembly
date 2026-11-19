process CRAMALIGN_GENCRAMCHUNKS {
    label "process_single"
    executor "local"

    input:
    // Native processes can't take path values as inputs
    tuple val(meta), val(cram), val(crai)
    val cram_bin_size

    output:
    tuple val(meta), val(cram), val(crai), val(chunkn), val(slices), emit: cram_slices
    tuple val("${task.process}"), val('cramchunks'), val('1.1.0'), emit: versions_cramchunks, topic: versions

    when:
    task.ext.when == null || task.ext.when

    exec:
    def n_slices = file(crai).countLines(decompress: true)
    def size     = cram_bin_size
    def n_bins   = Math.ceil(n_slices / size).toInteger()
    chunkn       = (0..<n_bins).collect()
    slices       = chunkn.collect { chunk ->
        def lower = chunk * size
        def upper = [lower + size, n_slices].min()

        return [ lower, upper - 1 ]
    }
}
