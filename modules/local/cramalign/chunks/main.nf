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
    def slices   = (0..n_bins).collect { chunk ->
        def lower = chunk == 0 ? 0 : (chunk * size) + 1
        def upper = chunk == n_bins ? n_slices : (chunk + 1) * size
        return [ lower, upper ]
    }
    def chunkn = (0..n_bins)
}
