/*
    Publishing functions for sanger-tol/genomeassembly
*/

//
// Get Shane's short-style name for read platform
//
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

//
// Convert internal stage names to output dirs
//
def stageToAssemblyDir(stage) {
    if (stage == "assembly") { return "raw" }
    if (stage == "purging") { return "purging" }
    if (stage == "polishing") { return "polishing" }
    if (stage == "scaffolding") { return "scaffolding" }
}
