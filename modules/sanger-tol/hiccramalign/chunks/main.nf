process HICCRAMALIGN_CHUNKS {
    label "process_single"
    executor "local"

    input:
    // Native processes can't take path values as inputs
    tuple val(meta), val(cram), val(crai)
    val cram_bin_size

    output:
    tuple val(meta), val(cram), val(crai), val(chunkn), val(slices), emit: cram_slices
    path("versions.yml")                                           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    exec:
    // Note: Manually bump version number when updating module
    def VERSION = "1.0.0"

    def n_slices = file(crai).countLines(decompress: true) - 1
    def size     = cram_bin_size
    def n_bins   = n_slices.intdiv(size)
    chunkn       = (0..n_bins)
    slices       = chunkn.collect { chunk ->
        def lower = chunk * size
        def upper = [lower + size - 1, n_slices - 1].min()

        return [ lower, upper ]
    }
    def versions_file = new File("${task.workDir}/versions.yml")
    versions_file.write(
        """
        HICCRAMALIGN_CHUNKS:
            hiccramalign_chunks: ${VERSION}
        """
    )
}
