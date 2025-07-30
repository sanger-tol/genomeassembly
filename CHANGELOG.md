# sanger-tol/genomeassembly: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v0.10.1dev - [date]

### `Added`

Overhaul pipeline input YAML structure (by @prototaxites):

- Validate the input YAML file with nf-schema
- Re-structure the input YAML a bit to be clearer
  - Move all species-level metadata to a "metadata" subsection
  - Move all input sequence data to a "sequencing_data" subsection
    - Each sequence data readset has an optional FastK input section (directory, basename, kmer size)
    - The long reads have an optional "coverage" input allowing pre-specification of the 1n coverage
  - Move the database stuff (busco + oatk) to a "databases" subsection

Overhaul kmer subworkfow (by @prototaxites):

- Pass any read set (currently long reads, maternal + paternal reads) without a pre-supplied FastK db through FastK
- For long reads, estimate coverage with Genomescope2 only if coverage is not supplied
- For maternal and paternal reads, build Yak DBs and make the hapmer dbs using merquryfk/hapmaker

Overhaul raw assembly subworkflow (by @prototaxites):

- Hifiasm is first run in --bin-only mode to generate .bin files of error corrected reads
- These .bin files are then used to populate all assemblies
- Add a --hifiasm_options parameter to tune hifiasm command line arguments
- Pipeline then uses channel logic to decide
  - run base assembly (always) - if the assembly is run without --primary, then hap1/hap2 are selected - otherwise pri/alt are used
  - if hic_phasing is enabled AND hic data is available - run hic_phased assembly
  - if trio_binning is enabled AND trio data is available - run trio_binned assembly
- Replace GFA_TO_FASTA with GAWK

Overhaul purging subworkflow (by @prototaxites):

- Rename PURGE_DUPS to PURGING
- Add parameter purging_assemblytypes which takes a comma-separated list of types of assembly ("primary", "hic_phased", "trio_binned") to define which assembly types to purge. (default "primary")
- Add parameter purging_purge_middle which if set removes the -e arg from the purgedups/getseqs command, allowing purging in the middle of contigs
- Add parameter purging_cutoffs which if set overrides the automatic setting of the purge_dups calcuts cutoff parameters. Must be a comma-separated string with three entries - e.g. "5,20,100".
- Remove GET_CALCUTS_PARAMS and calculate these args directly from the coverage meta value
- Move the re-concatentation of purged haplotigs and the alternate assembly inside the PURGING subworkflow

Overhaul polishing subworkflow (by @prototaxites):

- Subworkflow can take multiple assemblies for polishing
- Tidy up code to use modern Nextflow practices; avoid name collisions using meta objects for prefix setting
- Replace BED_CHUNKS module with GAWK
- Remove LONGRANGER_COVERAGE and just pull this value out the summary file with some Groovy
- Move the hap1/hap2 splitting inside the subworkflow and accomplish this with Seqkit grep + known regular expressions
- Update all nf-core modules used by the subworkflow

Overhaul Hi-C mapping subworkflow (by @prototaxites):

- Use new sanger-tol shared Hi-C mapping subworkflow
- This can map with both bwa-mem2 or minimap + implements CRAM chunking

Overhaul scaffolding subworkflow (by @prototaxites):

- Update YaHS to newest version
- TODO: Move Hi-C contact map generation to new subworkflow(s)

### `Fixed`

### `Dependencies`

### `Deprecated`

## [[0.10.0](https://github.com/sanger-tol/genomeassembly/releases/tag/0.10.0)] - Hideous Zippleback - [2024-04-16]

### Enhancements & fixes

- OATK module is added into the ORGANELLES subworkflow
- ORGANELLES subworkflow is now called once in the main workflow and runs MITOHIFI in read and assembly mode along with OATK
- ORGANELLES module is now tested in github CI
- NCBI API secret introduced to run MITOHIFI_FINDMITOREFERENCE module
- hifiasm haplotigs are not purged anymore
- Longranger container version is updated

### Software dependencies

Note, since the pipeline is using Nextflow DSL2, each process will be run with its own [Biocontainer](https://biocontainers.pro/#/registry). This means that on occasion it is entirely possible for the pipeline to be using different versions of the same tool. However, the overall software dependency changes compared to the last release have been listed below for reference.

| Dependency | Old version | New version |
| ---------- | ----------- | ----------- |
| mitohifi   | 3.0.0       | 3.1.1       |
| oatk       |             | 1.0         |

**NB:** Dependency has been **added** if just the new version information is present.

## [[0.9.0](https://github.com/sanger-tol/genomeassembly/releases/tag/0.9.0)] - Night Fury - [2023-12-15]

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
| ------------- | ---------------------- |
|               | --input                |
|               | --bed_chunks_polishing |
|               | --cool_bin             |
|               | --polishing_on         |
|               | --hifiasm_hic_on       |
|               | --organelles_on        |

> **NB:** Parameter has been **added** if just the new parameter information is present.

### Software dependencies

Note, since the pipeline is using Nextflow DSL2, each process will be run with its own [Biocontainer](https://biocontainers.pro/#/registry). This means that on occasion it is entirely possible for the pipeline to be using different versions of the same tool. However, the overall software dependency changes compared to the last release have been listed below for reference.

| Dependency      | Old version | New version                              |
| --------------- | ----------- | ---------------------------------------- |
| bcftools        |             | 1.17, 1.16                               |
| bed_chunks.sh   |             | 1.0                                      |
| bedtools        |             | 2.31.0                                   |
| busco           |             | 5.4.3                                    |
| bwamem2         |             | 2.2.1                                    |
| coverage.sh     |             | 1.0                                      |
| cooler          |             | 0.8.11                                   |
| cut             |             | 8.30                                     |
| fastk           |             | f18a4e6d2207539f7b84461daebc54530a9559b0 |
| freebayes       |             | 1.3.6                                    |
| gatk4           |             | 4.4.0.0                                  |
| genescope       |             | 380815c420f50171f9234a0fd1ff426b39829b91 |
| gfastats        |             | 1.3.5                                    |
| GNU Awk         |             | 5.1.0                                    |
| hifiasm         |             | 0.19.3-r572                              |
| juicer          |             | 1.1                                      |
| juicer tools    |             | 1.9.9                                    |
| longranger      |             | 2.2.2                                    |
| merquryfk       |             | 8ae344092df5dcaf83cfb7f90f662597a9b1fc61 |
| minimap2        |             | 2.24-r1122                               |
| mitohifi        |             | 3.0.0                                    |
| pigz            |             | 2.3.4                                    |
| pretextmap      |             | 0.1.9                                    |
| pretextsnapshot |             | 0.0.4                                    |
| purge_dups      |             | 1.2.5                                    |
| python          |             | 3.11.0, 3.8.1                            |
| R               |             | 4.2.0                                    |
| samtools        |             | 1.17, 1.16.1                             |
| seqtk           |             | 1.3-r106                                 |
| yahs            |             | 1.2a.1                                   |
| yaml            |             | 6.0                                      |

**NB:** Dependency has been **added** if just the new version information is present.
