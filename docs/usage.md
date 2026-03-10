# sanger-tol/genomeassembly: Usage

## Introduction

## Input files

The pipeline requires two input files to run: a genomic data file and an assembly specifications file. These are provided using the `--genomic_data` and `--assembly_specs` parameters respectively, and can be in either JSON or YAML format.

```
--genomic_data '[path to genomic data file]' --assembly_specs '[path to assembly specs file]'
```

### Genomic data input

The genomic data input file describes the sequencing datasets available for assembly. Each entry in the array describes a dataset, including the sample identifier and the sequencing platform. The following fields are required:

| Field      | Required | Description                                                                                                                           |
|------------|----------|---------------------------------------------------------------------------------------------------------------------------------------|
| `id`       | Yes      | Sample identifier for the dataset.                                                                                                    |
| `platform` | Yes      | Sequencing platform. One of: `pacbio_clr`, `pacbio_hifi`, `oxford_nanopore`, `illumina`, `illumina_hic`, `illumina_10x`               |
| `reads`    | Yes      | Array of file paths to sequencing reads in FASTA, FASTQ, or CRAM format. The required file format depends on the sequencing platform. |
| `fastk`    | No       | Pre-computed FastK database (see below)                                                                                               |

The filetype of the input reads depends on the sequencing platform:

| Platform          | Allowed File Extensions                      |
|-------------------|----------------------------------------------|
| `pacbio_hifi`     | `.fna`, `.fa`, `.fasta` (optionally gzipped) |
| `oxford_nanopore` | `.fq`, `.fastq` (optionally gzipped)         |
| `illumina_hic`    | `.cram`                                      |
| `illumina_10x`    | `.fq`, `.fastq` (optionally gzipped)         |
| `illumina`        | `.cram                                       |

It is also possible to provide pre-computed FastK databases if desired, to skip computation by the pipeline. To to this, the FastK input requires all the following sub-fields:

| Field       | Description                                                        |
|-------------|--------------------------------------------------------------------|
| `hist`      | Path to the FastK histogram `.hist` file                           |
| `ktab`      | Array of paths to FastK ktab files, including hidden files         |
| `kmer_size` | Integer kmer size used to generate the FastK database (minimum: 5) |

#### Example genomic data file

```json
[
  {
    "id": "sample1",
    "platform": "pacbio_hifi",
    "reads": [
      "/path/to/hifi_reads1.fq.gz",
      "/path/to/hifi_reads2.fq.gz"
    ]
  },
  {
    "id": "sample1",
    "platform": "illumina_hic",
    "reads": [
      "/path/to/hic_reads_R1.fq.gz",
      "/path/to/hic_reads_R2.fq.gz"
    ]
  },
  {
    "id": "sample1",
    "platform": "illumina_10x",
    "reads": [
      "/path/to/10x_reads_R1.fq.gz",
      "/path/to/10x_reads_R2.fq.gz"
    ],
    "fastk": {
      "hist": "/path/to/fastk/sample.hist",
      "ktab": [
        "/path/to/fastk/sample.ktab",
        "/path/to/fastk/.sample.ktab.1"
      ],
      "kmer_size": 31
    }
  },
  {
    "id": "maternal_sample",
    "platform": "illumina",
    "reads": [
      "/path/to/maternal_reads.fq.gz"
    ]
  },
  {
    "id": "paternal_sample",
    "platform": "illumina",
    "reads": [
      "/path/to/paternal_reads.fq.gz"
    ]
  }
]
```

### Assembly specification input

The assembly specifications file defines what assemblies to produce, including the required datasets,
options to enable or disable specific pipeline stages, and parameters to tune the various pipeline tools. It
can be in either YAML or JSON format. The following fields are mandatory for all specifications:

| Parameter            | Description                                                  |
|----------------------|--------------------------------------------------------------|
| `id`                 | Unique identifier for the assembly.                          |
| `assembler`          | Which assembler to use: `hifiasm`, `mitohifi`, or `oatk`     |
| `long_read_dataset`  | Dataset name from the genomic data file supplying long reads |
| `long_read_platform` | Platform for long reads: `pacbio_hifi` or `oxford_nanopore`  |

#### Hifiasm assembly options

These options can be used if the `assembler` field is set to "hifiasm".

| Parameter                       |            | Description                                                                                                                              |
|---------------------------------|------------|------------------------------------------------------------------------------------------------------------------------------------------|
| `ultralong_dataset`             | -          | Sample identifier from which to find ultra-long reads for Hifiasm assembly. Must have an `oxford_nanopore` platform entry.               |
| `hic_dataset`                   | -          | Sample identifier from which to find Hi-C reads for phased Hifiasm assembly and scaffolding. Must have an `illumina_hic` platform entry. |
| `polishing_dataset`             | -          | Sample identifier from which to find 10X reads for polishing. Must have an `illumina_10x` platform entry.                                |
| `maternal_dataset`              | -          | Sample identifier from which to find maternal reads for trio assembly.                                                                   |
| `maternal_platform`             | -          | Platform for maternal reads: `illumina`, `illumina_10x`, `pacbio_hifi`, or `oxford_nanopore`                                             |
| `paternal_dataset`              | -          | Sample identifier from which to find maternal reads for trio assembly.                                                                   |
| `paternal_platform`             | -          | Platform for paternal reads: `illumina`, `illumina_10x`, `pacbio_hifi`, or `oxford_nanopore`                                             |
| `long_read_1n_coverage`         | -          | Haploid/1n coverage of the target genome in the long read dataset                                                                        |
| `phased_assembly`               | `false`    | Produce a phased Hifiasm assembly with Hi-C data. Requires `hic_dataset`.                                                                |
| `trio_assembly`                 | `false`    | Produce a trio-binned Hifiasm assembly. Requires `maternal_dataset` and `paternal_dataset`.                                              |
| `purge`                         | `false`    | Purge retained haplotypic duplications using the purge_dups pipeline.                                                                    |
| `polish`                        | `false`    | Polish the assembly using Illumina 10X data                                                                                              |
| `scaffold`                      | `true`     | Map Hi-C reads and scaffold the assembly using YaHS.                                                                                     |
| `find_mito`                     | `true`     | Enable mitochondrial genome assembly/search with Mitohifi.                                                                               |
| `find_plastid`                  | `false`    | Enable plastid genome assembly/search with Mitohifi.                                                                                     |
| `hifiasm_bin_arguments`         | -          | Additional arguments for Hifiasm overlap graph generation (e.g., error correction options)                                               |
| `hifiasm_arguments`             | -          | Additional arguments for Hifiasm to produce an assembly.                                                                                 |
| `purging_cutoffs`               | -          | Comma-separated coverage cutoffs for purging (e.g., `"5,20,100"`).                                                                       |
| `purge_middle`                  | `false`    | Purge haplotypic duplications from within contigs, not just at the ends.                                                                 |
| `yahs_arguments`                | -          | Additional arguments for YaHS                                                                                                            |
| `busco_lineage`                 | `auto_euk` | BUSCO lineage for completeness assessment (e.g., `metazoa_odb12`)                                                                        |
| `mitohifi_reference_species`    | -          | Binomial name of taxon for reference mitochondrial genome. Required if `find_mito` or `find_plastid` are set.                            |
| `mitohifi_mito_genetic_code`    | -          | Mitochondrial genetic code for gene prediction. Required if `find_mito` or `find_plastid` are set.                                       |
| `mitohifi_plastid_genetic_code` | 11         | Plastid genetic code for gene prediction. Required if `find_mito` or `find_plastid` are set.                                             |
| `mitohifi_arguments`            | -          | Extra arguments for MitoHiFi in reads mode                                                                                               |

#### MitoHiFi options

The following options can be used when `assembler` is "mitohifi". Currently, this requires that `long_read_platform` is "pacbio_hifi".

| Parameter                    | Default | Description                                                           |
|------------------------------|---------|-----------------------------------------------------------------------|
| `mitohifi_reference_species` | -       | Binomial name of taxon for reference mitochondrial genome. (Required) |
| `mitohifi_mito_genetic_code` | -       | Mitochondrial genetic code for gene prediction. (Required)            |
| `mitohifi_arguments`         | -       | Extra arguments for MitoHiFi in reads mode                            |

#### Oatk options

The following options can be used when `assembler` is "mitohifi". Currently, this requires that `long_read_platform` is "pacbio_hifi". One or both of `oatk_mito_hmm` or `oatk_plastid_hmm` is required.

| Parameter              | Default | Description                                                                    |
|------------------------|---------|--------------------------------------------------------------------------------|
| `oatk_kmer_size`       | `1000`  | Kmer size for oatk                                                             |
| `oatk_coverage_cutoff` | -       | Coverage cutoff for oatk. Auto-calculated as `(5 x coverage)` if not provided. |
| `oatk_arguments`       | -       | Additional arguments for oatk (excluding -k and -c)                            |
| `oatk_mito_hmm`        | -       | Path to oatk mitochondrial HMM file (`.fam` format)                            |
| `oatk_plastid_hmm`     | -       | Path to oatk plastid HMM file (`.fam` format)                                  |

#### Example assembly specifications file

```json
- id: sample1_assembly_hifiasm
  assembler: hifiasm
  long_read_dataset: sample1
  long_read_platform: pacbio_hifi
  hic_dataset: sample1
  scaffold: true
  phased_assembly: true
  purge: true
  find_mito: true
  mitohifi_reference_species: Bellis perennis
  mitohifi_mito_genetic_code: 1
  busco_lineage: eudicots_odb10
- id: sample1_assembly_trio
  assembler: hifiasm
  long_read_dataset: sample1
  long_read_platform: pacbio_hifi
  hic_dataset: sample1
  trio_assembly: true
  maternal_dataset: maternal_sample
  maternal_platform: illumina
  paternal_dataset: paternal_sample
  paternal_platform: illumina
  scaffold: true
  find_mito: false
- id: sample1_mito
  assembler: oatk
  long_read_dataset: sample1
  long_read_platform: pacbio_hifi
  oatk_mito_hmm: /path/to/mito.fam
```

## Additional setup procedures

### Longranger (polishing)

Longranger is a proprietary software product from 10X Genomics.
Its terms and conditions state that we _cannot_ redistribute the copy we use in the Tree of Life department.

To run the polishing subroutine of this pipeline, you will have to install longranger yourself.

Go to <https://support.10xgenomics.com/genome-exome/software/downloads/latest>,
read their End User Software License Agreement, and you'll be able to download the software if you accept it.

To make a Docker (or Singularity) container out of it, use the following Dockerfile.

```Dockerfile
FROM ubuntu:22.04
LABEL org.opencontainers.image.licenses="10x Genomics End User Software License Agreement - https://support.10xgenomics.com/genome-exome/software/downloads/latest"
ARG DEST=/opt
ADD ./longranger-2.2.2.tar.gz $DEST
RUN ln -s $DEST/longranger-2.2.2/longranger /usr/local/bin/
```

Then, to use the container in the pipeline, pass the path or name of the container to the
`--polishing_longranger_container_path` parameter when running the pipeline.

### NCBI API Key

Running Mitohifi for organelle assembly requires access to the NCBI API. Although this
is possible without configuration, it is possible to specify an NCBI API key, if you have
one. To do this, configure the [Nextflow secret](https://www.nextflow.io/docs/latest/secrets.html)
NCBI_API_KEY as follows:

```bash
  nextflow secrets set NCBI_API_KEY '[API key]'
```

## Running the pipeline

The typical command for running the pipeline is as follows:

```console
nextflow run sanger-tol/genomeassembly --genomic_data assets/genomic_data.json --assembly_specs assets/assembly_specs.json --outdir <OUTDIR> -profile docker
```

This will launch the pipeline with the `docker` configuration profile. See below for more information about profiles.

Note that the pipeline will create the following files in your working directory:

```bash
work                # Directory containing the nextflow working files
<OUTDIR>            # Finished results in specified location (defined with --outdir)
.nextflow_log       # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

If you wish to repeatedly use the same parameters for multiple runs, rather than specifying each flag in the command, you can specify these in a params file.

Pipeline settings can be provided in a `yaml` or `json` file via `-params-file <file>`.

> [!WARNING]
> Do not use `-c <file>` to specify parameters as this will result in errors. Custom config files specified with `-c` must only be used for [tuning process resource specifications](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources), other infrastructural tweaks (such as output directories), or module arguments (args).

The above pipeline run specified with a params file in yaml format:

```bash
nextflow run sanger-tol/genomeassembly -profile docker -params-file params.yaml
```

with:

```yaml title="params.yaml"
genomic_data: './genomic_data.json'
assembly_specs: './assembly_specs.json'
outdir: './results/'
<...>
```

You can also generate such `YAML`/`JSON` files via [sanger-tol/launch](https://pipelines.tol.sanger.ac.uk/launch).

### Updating the pipeline

When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. When running the pipeline after this, it will always use the cached version if available - even if the pipeline has been updated since. To make sure that you're running the latest version of the pipeline, make sure that you regularly update the cached version of the pipeline:

```console
nextflow pull sanger-tol/genomeassembly
```

### Reproducibility

It is a good idea to specify the pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.

First, go to the [sanger-tol/genomeassembly releases page](https://github.com/sanger-tol/genomeassembly/releases) and find the latest pipeline version - numeric only (eg. `1.3.1`). Then specify this when running the pipeline with `-r` (one hyphen) - eg. `-r 1.3.1`. Of course, you can switch to another version by changing the number after the `-r` flag.

This version number will be logged in reports when you run the pipeline, so that you'll know what you used when you look back in the future.

To further assist in reproducibility, you can use share and reuse [parameter files](#running-the-pipeline) to repeat pipeline runs with the same settings without having to write out a command with every single parameter.

> [!TIP]
> If you wish to share such profile (such as upload as supplementary material for academic publications), make sure to NOT include cluster specific paths to files, nor institutional specific profiles.

## Core Nextflow arguments

> [!NOTE]
> These options are part of Nextflow and use a _single_ hyphen (pipeline parameters use a double-hyphen)

### `-profile`

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments.

Several generic profiles are bundled with the pipeline which instruct the pipeline to use software packaged using different methods (Docker, Singularity, Podman, Shifter, Charliecloud, Apptainer, Conda) - see below.

> [!IMPORTANT]
> We highly recommend the use of Docker or Singularity containers for full pipeline reproducibility, however when this is not possible, Conda is also supported.

The pipeline also dynamically loads configurations from [https://github.com/nf-core/configs](https://github.com/nf-core/configs) when it runs, making multiple config profiles for various institutional clusters available at run time. For more information and to check if your system is supported, please see the [nf-core/configs documentation](https://github.com/nf-core/configs#documentation).

Note that multiple profiles can be loaded, for example: `-profile test,docker` - the order of arguments is important!
They are loaded in sequence, so later profiles can overwrite earlier profiles.

If `-profile` is not specified, the pipeline will run locally and expect all software to be installed and available on the `PATH`. This is _not_ recommended, since it can lead to different results on different machines dependent on the computer environment.

- `test`
  - A profile with a complete configuration for automated testing
  - Includes links to test data so needs no other parameters
- `docker`
  - A generic configuration profile to be used with [Docker](https://docker.com/)
- `singularity`
  - A generic configuration profile to be used with [Singularity](https://sylabs.io/docs/)
- `podman`
  - A generic configuration profile to be used with [Podman](https://podman.io/)
- `shifter`
  - A generic configuration profile to be used with [Shifter](https://nersc.gitlab.io/development/shifter/how-to-use/)
- `charliecloud`
  - A generic configuration profile to be used with [Charliecloud](https://charliecloud.io/)
- `apptainer`
  - A generic configuration profile to be used with [Apptainer](https://apptainer.org/)
- `wave`
  - A generic configuration profile to enable [Wave](https://seqera.io/wave/) containers. Use together with one of the above (requires Nextflow ` 24.03.0-edge` or later).
- `conda`
  - A generic configuration profile to be used with [Conda](https://conda.io/docs/). Please only use Conda as a last resort i.e. when it's not possible to run the pipeline with Docker, Singularity, Podman, Shifter, Charliecloud, or Apptainer.

### `-resume`

Specify this when restarting a pipeline. Nextflow will use cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously. For input to be considered the same, not only the names must be identical but the files' contents as well. For more info about this parameter, see [this blog post](https://www.nextflow.io/blog/2019/demystifying-nextflow-resume.html).

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

### `-c`

Specify the path to a specific config file (this is a core Nextflow command). See the [nf-core website documentation](https://nf-co.re/usage/configuration) for more information.

## Custom configuration

### Resource requests

Whilst the default requirements set within the pipeline will hopefully work for most people and with most input data, you may find that you want to customise the compute resources that the pipeline requests. Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the pipeline steps, if the job exits with any of the error codes specified [here](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/conf/base.config#L18) it will automatically be resubmitted with higher resources request (2 x original, then 3 x original). If it still fails after the third attempt then the pipeline execution is stopped.

To change the resource requests, please see the [max resources](https://nf-co.re/docs/usage/configuration#max-resources) and [tuning workflow resources](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources) section of the nf-core website.

### Custom Containers

In some cases, you may wish to change the container or conda environment used by a pipeline steps for a particular tool. By default, nf-core pipelines use containers and software from the [biocontainers](https://biocontainers.pro/) or [bioconda](https://bioconda.github.io/) projects. However, in some cases the pipeline specified version maybe out of date.

To use a different container from the default container or conda environment specified in a pipeline, please see the [updating tool versions](https://nf-co.re/docs/usage/configuration#updating-tool-versions) section of the nf-core website.

### Custom Tool Arguments

A pipeline might not always support every possible argument or option of a particular tool used in pipeline. Fortunately, nf-core pipelines provide some freedom to users to insert additional parameters that the pipeline does not include by default.

To learn how to provide additional arguments to a particular tool of the pipeline, please see the [customising tool arguments](https://nf-co.re/docs/usage/configuration#customising-tool-arguments) section of the nf-core website.

### nf-core/configs

In most cases, you will only need to create a custom config as a one-off but if you and others within your organisation are likely to be running nf-core pipelines regularly and need to use the same settings regularly it may be a good idea to request that your custom config file is uploaded to the `nf-core/configs` git repository. Before you do this please can you test that the config file works with your pipeline of choice using the `-c` parameter. You can then create a pull request to the `nf-core/configs` repository with the addition of your config file, associated documentation file (see examples in [`nf-core/configs/docs`](https://github.com/nf-core/configs/tree/master/docs)), and amending [`nfcore_custom.config`](https://github.com/nf-core/configs/blob/master/nfcore_custom.config) to include your custom profile.

See the main [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for more information about creating your own configuration files.

If you have any questions or issues please send us a message on [Slack](https://nf-co.re/join/slack) on the [`#configs` channel](https://nfcore.slack.com/channels/configs).

## Running in the background

Nextflow handles job submissions and supervises the running jobs. The Nextflow process must run until the pipeline is finished.

The Nextflow `-bg` flag launches Nextflow in the background, detached from your terminal so that the workflow does not stop if you log out of your session. The logs are saved to a file.

Alternatively, you can use `screen` / `tmux` or similar tool to create a detached session which you can log back into at a later time.
Some HPC setups also allow you to run nextflow within a cluster job submitted your job scheduler (from where it submits more jobs).

## Nextflow memory requirements

In some cases, the Nextflow Java virtual machines can start to request a large amount of memory.
We recommend adding the following line to your environment to limit this (typically in `~/.bashrc` or `~./bash_profile`):

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```
