/*
   sanger-tol/genomeassembly input validation functions
*/

// Validate that all the reads in a dataset have the same file extension, and that
// the file extensions are correct for the sequencing platform.
def validateReadFiles(meta, reads) {
    // Check that all files have the same extension
    def extensions = reads.collect { file ->
        if(file.getExtension() == "gz") {
            return file.getName().tokenize(".").takeRight(2).join('.')
        } else {
            return file.getExtension()
        }
    }.unique()

    if(extensions.size() > 1) {
        error("Dataset validation error [${meta.id}:${meta.platform}]: Not all files are of the same type!")
    }

    def platform_ext_map = [
        ["pacbio_hifi", /^fn?(asta)?(\.gz)?$/],
        ["oxford_nanopore", /^f(ast)?q(\.gz)?$/],
        ["illumina_hic", /^cram$/],
        ["illumina_10x", /^f(ast)?q(\.gz)?$/],
        ["illumina", /^cram$/],
    ].find { platform, _regex -> platform == meta.platform }

    // Validate that the correct extension is present for each platform
    if(!(extensions[0] =~ platform_ext_map[1])) {
        error("Dataset validation error [${meta.id}:${meta.platform}]: File extension ${extensions[0]} does not match the expected input ${platform_ext_map[1]}.")
    }
}

// For each dataset in an assembly specification, check that it exists in the data
// specification.
def checkDataExists(spec, datasets) {
    def platform_key = [
        [name: 'long_read_dataset', platform: spec.long_read_platform],
        [name: 'hic_dataset', platform: spec.hic_platform],
        [name: 'polishing_dataset', platform: spec.polishing_platform],
        [name: 'maternal_dataset', platform: spec.maternal_platform],
        [name: 'paternal_dataset', platform: spec.paternal_platform]
    ].findAll { dataKey -> spec[dataKey.name] }

    platform_key.each { platform ->
        if(spec[platform.name] && !datasets.find { data -> data.id == spec[platform.name] && data.platform == platform.platform } ) {
            error("Assembly specification error [${spec.id}]: There is no ${platform.platform} dataset called '${spec[platform.name]}' does not exist!")
        }
    }
}
