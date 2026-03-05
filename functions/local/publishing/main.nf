/**
 * Output Publishing Functions
 *
 * Transforms internal naming conventions to output-ready formats: platform names
 * to abbreviations and assembly stages to output directory structures.
 */

/**
 * Converts sequencing platform names to short abbreviations for file naming.
 *
 * @param platform Full platform name (pacbio_hifi, oxford_nanopore, illumina_10x, illumina_hic, illumina)
 *
 * @return Short platform abbreviation (ccs, ont, 10x, illumina) or empty string if unknown
 */
def getPlatformShortName(platform) {
    def shortname = ""
    if(platform == "pacbio_hifi") {
        shortname = "ccs"
    } else if (platform == "oxford_nanopore") {
        shortname = "ont"
    } else if (platform == "illumina_10x"){
        shortname = "10x"
    } else if (platform == "illumina_hic") {
        shortname == "illumina"
    } else if (platform == "illumina") {
        shortname == "illumina"
    }

    return shortname
}

/**
 * Maps internal assembly stage names to their corresponding output directory names.
 *
 * @param stage Internal stage name (assembly, purging, polishing, scaffolding)
 *
 * @return Output directory name (raw, purging, polishing, scaffolding) or null if unknown
 */
def specToAssemblyDir(spec) {
    def stage = spec.stage

    if(spec.assembler == "hifiasm") {
        if (stage == "assembly") { return "raw" }
        if (stage == "purging") { return "purging" }
        if (stage == "polishing") { return "polishing" }
        if (stage == "scaffolding") { return "scaffolding" }
        if (stage == "mitohifi_mito") { return "mito" }
        if (stage == "mitohifi_plastid") { return "plastid" }
    } else if (spec.assembler == "mitohifi") {
        return ""
    } else if (spec.assembler == "oatk") {
        return ""
    }
}
