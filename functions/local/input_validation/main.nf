/**
 * Input validation for sequencing datasets and assembly specifications.
 *
 * Functions validate:
 * - Read file consistency (extensions match platform expectations)
 * - Dataset existence and availability
 * - HMMER3 database file integrity
 */

/**
 * Validates that all read files in a dataset have consistent file extensions and that
 * those extensions match the expected format for the sequencing platform.
 *
 * This function performs two key validations:
 * 1. Ensures all files in the dataset share the same extension (e.g., all .fastq.gz or all .cram)
 * 2. Verifies the extension matches the platform's expected format (e.g., PacBio HiFi uses .fasta)
 *
 * Supported platforms and their expected extensions:
 * - pacbio_hifi: .fna, .fa, .fasta (optionally .gz compressed)
 * - oxford_nanopore: .fq, .fastq (optionally .gz compressed)
 * - illumina_hic: .cram (uncompressed only)
 * - illumina_10x: .fq, .fastq (optionally .gz compressed)
 * - illumina: .cram (uncompressed only)
 *
 * @param meta Metadata map containing dataset identifiers:
 *             [
 *               id: <dataset identifier>
 *               platform: <sequencing platform name>
 *             ]
 * @param reads List of file objects representing the sequencing read files to validate
 *
 * @throws Error if files have inconsistent extensions or if extensions don't match
 *         the platform requirements. Error messages include the dataset id and platform
 *         for debugging purposes.
 */
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
        ["pacbio_hifi", /^fn?(ast)?a(\.gz)?$/],
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

/**
 * Validates that all datasets referenced in an assembly specification exist in the
 * available datasets list.
 *
 * This function checks that each dataset referenced by the assembly specification
 * (identified by dataset name and platform) can be found in the provided datasets list.
 *
 * Only datasets that are specified (non-null) in the assembly spec are validated.
 *
 * @param spec The assembly specification containing dataset references in the form
 *             of <type>_dataset (dataset name) and <type>_platform (platform name) pairs.
 *             Must also contain an 'id' field for error reporting.
 * @param datasets List of available dataset objects, each containing at minimum:
 *                 [
 *                   id: <dataset identifier>
 *                   platform: <sequencing platform name>
 *                 ]
 *
 * @throws Error if any referenced dataset cannot be found in the datasets list.
 *         Error messages include the assembly id, platform, and dataset name for debugging.
 */
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

/**
 * Collects all required HMMER3 database files associated with a base HMM profile.
 *
 * Returns the base HMM file plus its four required index files (.h3f, .h3i, .h3m, .h3p).
 *
 * @param hmmPath Path to the base HMM profile file (without index extensions)
 *
 * @return List of file objects: [base HMM file, .h3f, .h3i, .h3m, .h3p]
 *
 * @throws FileNotFoundException if any required file cannot be found
 */
def createHmmFilesList(hmmPath) {
    def hmm_extensions = [".h3f", ".h3i", ".h3m", ".h3p"]

    return [file(hmmPath, checkIfExists: true)] + hmm_extensions.collect { ext -> file(hmmPath + ext, checkIfExists: true) }
}
