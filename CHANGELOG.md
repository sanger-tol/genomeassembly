# sanger-tol/genomeassembly: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [[1.0.0](https://github.com/sanger-tol/genomeassembly/releases/tag/1.0.0)] - Night Fury - [2023-12-15]

Initial release of sanger-tol/genomeassembly, created with the [nf-core](https://nf-co.re/) template.

### Enhancements & fixes

- Created with nf-core/tools template v2.8.0.
- Implements the basic genome assembly pipeline used in ToL
- Subworkflow to generate Genomescope model
- Subworkflow to generate a raw genome assembly
- Subworkflow to purge haplotigs from primary
- Subworkflow for polishing with 10X data
- Subworkflow for HiC mapping
- Subworkflow for scaffolding
- Subworkflow to run organelles assembly
- Subworkflow to estimate statistics using gfastats, BUSCO, MerquryFK.

### Parameters

| Old parameter | New parameter          |
| ------------- | -----------------------|
|               | --input                |
|               | --bed_chunks_polishing |
|               | --cool_bin             |
|               | --polishing_on         |
|               | --hifiasm_hic_on       |
|               | --organelles_on        |

> **NB:** Parameter has been **added** if just the new parameter information is present. 

### Software dependencies

Note, since the pipeline is using Nextflow DSL2, each process will be run with its own [Biocontainer](https://biocontainers.pro/#/registry). This means that on occasion it is entirely possible for the pipeline to be using different versions of the same tool. However, the overall software dependency changes compared to the last release have been listed below for reference.

| Dependency        | Old version | New version                              |
| ----------------- | ----------- | ---------------------------------------- |
| bcftools          |             | 1.17, 1.16                               |
| bed_chunks.sh     |             | 1.0                                      |
| bedtools          |             | 2.31.0                                   |
| busco             |             | 5.4.3                                    |
| bwamem2           |             | 2.2.1                                    |
| coverage.sh       |             | 1.0                                      |
| cooler            |             | 0.8.11                                   |
| cut               |             | 8.30                                     |
| fastk             |             | f18a4e6d2207539f7b84461daebc54530a9559b0 |
| freebayes         |             | 1.3.6                                    |
| gatk4             |             | 4.4.0.0                                  |
| gfastats          |             | 1.3.5                                    |
| GNU Awk           |             | 5.1.0                                    |
| hifiasm           |             | 0.19.3-r572                              |
| juicer            |             | 1.1                                      |
| juicer tools      |             | 1.9.9                                    |
| longranger        |             | 2.2.2                                    |
| merquryfk         |             | 8ae344092df5dcaf83cfb7f90f662597a9b1fc61 |
| minimap2          |             | 2.24-r1122                               |
| mitohifi          |             | 3.0.0                                    |
| pigz              |             | 2.3.4                                    |
| pretextmap        |             | 0.1.9                                    |
| pretextsnapshot   |             | 0.0.4                                    |
| purge_dups        |             | 1.2.5                                    |
| python            |             | 3.11.0, 3.8.1                            |
| R                 |             | 4.2.0                                    |
| samtools          |             | 1.17, 1.16.1                             |
| seqtk             |             | 1.3-r106                                 |
| yahs              |             | 1.2a.1                                   |
| yaml              |             | 6.0                                      |

**NB:** Dependency has been **added** if just the new version information is present. 
