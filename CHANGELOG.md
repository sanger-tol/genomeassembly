# sanger-tol/genomeassembly: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v0.50.0 - TBD - [2025-12-10]

This is an interrim development release, tagging all the major changes to the pipeline since the last release, and
is likely the final release before the release of version 1.0.

In this release, the main codebase of the pipeline has been heavily overhauled to be more stable, and to meet newer
Nextflow standards, while trying to maintain general compatibility with the interface of the previous version, where possible.
However, note that the input YAML format has been changed heavily, so old YAML input files will need to be re-written.

The full 1.0 release will encompass further changes, particularly another overhaul of the input schema such that the data
input schema is separated from the assembly input schema, allowing multiple assemblies to be parameterised simultaneously
using specific data resources.

The full changelog is as follows:

### `Added`

General:

- Genericise file output prefixies to `asm` or `asm.{hap}` (by @prototaxites)
- Update nf-core template to 3.5.1 (by @prototaxites).

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

- Use sanger-tol/fasta_purge_retained_haplotype subworkflow
  - Add parameter purging_reads_chunk_size which defines how many reads to map in a single chunked job
  - Add parameter purging_assemblytypes which takes a comma-separated list of types of assembly ("primary", "hic_phased", "trio_binned") to define which assembly types to purge. (default "primary")
  - Add parameter purging_purge_middle which if set removes the -e arg from the purgedups/getseqs command, allowing purging in the middle of contigs
  - Add parameter purging_cutoffs which if set overrides the automatic setting of the purge_dups calcuts cutoff parameters. Must be a comma-separated string with three entries - e.g. "5,20,100".
  - Remove GET_CALCUTS_PARAMS and calculate these args directly from the coverage meta value

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
- only runs if --enable_scaffolding set
- use nf-core bam_stats_samtools subworkflow for mapping statistics

Overhaul scaffolding subworkflow (by @prototaxites):

- use new sanger-tol scaffolding subworkflow
  - Update YaHS to newest version
  - new module YAHS_MAKEPAIRSFILE which creates a pairs file of reads mapping to scaffolds - replaces PREPARE_PRETEXTMAP_INPUT
  - Update PretextMap + PretextSnapshot to latest versions
  - Update Cooler modules, overhauling their inputs
  - Use new JUICERTOOLS_PRE nf-core module, so we don't have to ship the jar any more
    - Supports multi-threading
    - This module has `errorStragegy = "ignore"` set to skip failures.

Overhaul genome statistics subworkflow (by @prototaxites):

- Use new sanger-tol shared genome statistics subworkflow
- Update to BUSCO 6.0

Overhaul organelle assembly subworkflow (by @prototaxites)

- Added parameters to control behaviour of oatk and mitohifi
  - sets `errorStragegy = "ignore"` for mitohifi/findmitoreference and mitohifi/mitohifi itself
  - `--oatk_kmer_size` and `--oatk_coverage`
  - `--mitohifi_min_ref_len`, `--mitohifi_reads_args` and `--mitohifi_contigs_args`

### `Fixed`

### `Dependencies`

| Module                     | Tool         | Old version                              | New version    |
| -------------------------- | ------------ | ---------------------------------------- | -------------- |
| bcftools/concat            | bcftools     | 1.16                                     | 1.21           |
| bcftools/consensus         | bcftools     | 1.16                                     | 1.21           |
| bcftools/index             | bcftools     | 1.16                                     | 1.21           |
| bcftools/norm              | bcftools     | 1.16                                     | 1.21           |
| bcftools/sort              | bcftools     | 1.16                                     | 1.21           |
| bcftools/view              | bcftools     | 1.16                                     | 1.21           |
| busco/busco                | busco        | 5.4.3                                    | 6.0.0          |
| bwamem2/index              | bwa-mem2     | 2.2.1                                    | 2.3            |
| cooler/cload               | cooler       | 0.8.11                                   | 0.10.4         |
| cooler/zoomify             | cooler       | 0.8.11                                   | 0.10.4         |
| fastk/fastk                | fastk        | f18a4e6d2207539f7b84461daebc54530a9559b0 | 1.1.0          |
| fastk/histex               | fastk        | f18a4e6d2207539f7b84461daebc54530a9559b0 | 1.1.0          |
| gatk4/mergevcfs            | gatk4        | 4.4.0.0                                  | 4.6.1.0        |
| gatk4/mergevcfs            | gcnvkernel   | <unspecified>                            | 0.9            |
| gawk                       | gawk         | -                                        | 5.3.0          |
| genescopefk                | genescopefk  | 380815c420f50171f9234a0fd1ff426b39829b91 | -              |
| genomescope2               | genomescope2 | -                                        | 2.0            |
| gfastats                   | gfastats     | 1.3.5                                    | 1.3.10         |
| juicertools/pre            | juicertools  | 1.9.9                                    | 2.20.00        |
| merquryfk/hapmaker         | merquryfk    | -                                        | 1.1.1          |
| merquryfk/merquryfk        | merquryfk    | 8ae344092df5dcaf83cfb7f90f662597a9b1fc61 | 1.1.1          |
| minimap2/align             | minimap2     | 2.24                                     | 2.30           |
| minimap2/align             | samtools     | 1.14                                     | 1.21           |
| minimap2/index             | minimap2     | 2.24                                     | 2.30           |
| mitohifi/findmitoreference | mitohifi     | 3.2.0                                    | 3.2.3          |
| mitohifi/mitohifi          | mitohifi     | 3.2.0                                    | 3.2.3          |
| pretextmap                 | pretextmap   | 0.1.8                                    | 0.1.9          |
| pretextmap                 | samtools     | 1.14                                     | 1.17           |
| purgedups/histplot         | purgedups    | -                                        | 1.2.6          |
| samtools/faidx             | samtools     | 1.17                                     | 1.22.1         |
| samtools/flagstat          | samtools     | 1.17                                     | 1.22.1         |
| samtools/idxstats          | samtools     | 1.17                                     | 1.22.1         |
| samtools/index             | samtools     | 1.17                                     | 1.22.1         |
| samtools/merge             | samtools     | 1.17                                     | 1.22.1         |
| samtools/splitheader       | samtools     | 1.17                                     | 1.22.1         |
| samtools/stats             | samtools     | 1.17                                     | 1.22.1         |
| seqkit/grep                | seqkit       | -                                        | 2.9.0          |
| tabix/bgzip                | tabix        | -                                        | 1.21           |
| yahs                       | yahs         | 1.2a.2                                   | 1.2.2          |
| yak/count                  | yak          | -                                        | 0.1            |
| asmstats                   | asmstats     | -                                        | 1.0            |
| bedtools/bamtobedsort      | bedtools     | -                                        | 2.31.1         |
| cramalign/bwamem2alignhic  | bwa-mem2     | -                                        | 2.2.1          |
| cramalign/bwamem2alignhic  | samtools     | -                                        | 1.22.1         |
| cramalign/minimap2alignhic | minimap2     | -                                        | 2.30           |
| cramalign/minimap2alignhic | samtools     | -                                        | 1.22.1         |
| fastxalign/minimap2align   | pyfastx      | -                                        | 2.2.0          |
| fastxalign/minimap2align   | minimap2     | -                                        | 2.30           |
| fastxalign/pyfastxindex    | pyfastx      | -                                        | 2.2.0          |
| hifiasm                    | hifiasm      | 0.18.5                                   | 0.25.0         |
| longranger/align           | longranger   | 2.2.2                                    | <user defined> |
| longranger/mkref           | longranger   | 2.2.2                                    | <user defined> |
| samtools/mergedup          | samtools     | -                                        | 1.22.1         |
| yahs/makepairsfile         | yahs         | -                                        | 1.2.2          |

### `Deprecated`

`--hifiasm_hic_on` replaced with `--enable_hic_phasing`
`--organelles_on` replaced with `--enable_organelle_assembly`
`--polishing_on` replaced with `--polishing_assemblytypes`
`--bed_chunks_polishing` replaced with `--polishing_n_chunks_bed`

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
