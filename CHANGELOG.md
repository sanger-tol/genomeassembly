# sanger-tol/genomeassembly: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [[0.9.1](https://github.com/sanger-tol/genomeassembly/releases/tag/0.9.1)] - Hideous Zippleback - [2024-04-15]

Initial release of sanger-tol/genomeassembly, created with the [nf-core](https://nf-co.re/) template.

### Enhancements & fixes

- OATK module is added into the ORGANELLES subworkflow 
- ORGANELLES subworkflow is now called once in the main workflow and runs MITOHIFI in read and assembly mode along with OATK
- ORGANELLES module is now tested in github CI
- NCBI API secret introduced to run MITOHIFI_FINDMITOREFERENCE module
- hifiasm haplotigs are not purged anymore
- Longranger container version is updated

### Software dependencies

Note, since the pipeline is using Nextflow DSL2, each process will be run with its own [Biocontainer](https://biocontainers.pro/#/registry). This means that on occasion it is entirely possible for the pipeline to be using different versions of the same tool. However, the overall software dependency changes compared to the last release have been listed below for reference.

| Dependency      | Old version | New version                              |
| --------------- | ----------- | ---------------------------------------- |
| mitohifi        | 3.0.0       | 3.1.1                                    |
| oatk            |             | 1.0                                      |


**NB:** Dependency has been **added** if just the new version information is present.
