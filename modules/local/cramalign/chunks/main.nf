process CRAM_CHUNKS {
    executor "local"

    input:
    tuple val(meta), val(cram), val(crai)
    val cram_bin_size

    output:
    tuple val(meta), val(cram), val(crai), val(chunkn), val(slices), emit: cram_slices

    exec:
    def n_slices = file(crai).countLines(decompress: true) - 1
    def size     = cram_bin_size
    def n_bins   = n_slices.intdiv(size)
    def chunkn = (0..n_bins)
    def slices   = chunkn.collect { chunk ->
        def lower = chunk * size
        def upper = [lower + size - 1, n_slices - 1].min()

        return [ lower, upper ]
    }
}
