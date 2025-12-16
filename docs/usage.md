# sanger-tol/genomeassembly: Usage

## Introduction

## Samplesheet input

You will need to create a samplesheet with information about the data you would like to assemble before running the pipeline. Use this parameter to specify its location. It has to be a YAML file, following the specification in the [Full samplesheet](#full-samplesheet) section below.

```
--input '[path to samplesheet file]'
```

### Full samplesheet

The input data for the pipeline is described in [YAML format](https://yaml.org/). It is a human readable file which contains information
about location paths for the raw data (HiFi, 10X, HiC) used for the genome assembly. It can also contain meta information such as HiC restriction motifs,
BUSCO lineage, mitochondrial code etc.

Briefly, the YAML describes three main subsections. The first is `metadata`, which described metadata about the target species, including a sample ID, the species binomial, and the mitochondrial genetic code. The second is `sequencing_data`, which describes the locations of the sequencing data to use for assembly. This has four sections (`long_reads`, `hic`, `illumina_10x`, `trio`), of which the only required section is `long_reads`. Any specified section **must** have a `reads` section, and the `long_reads` and `hic` sections also have a required `platform` and `hic_motif` section, respectively. All sections can optionally have a FastK database specified in the `fastk` section. Finally, there is a `databases` section which specifies which BUSCO and `oatk` databases to use.

An list of all the possible fields is given below, with their contents described as comments. Optional fields are tagged with `(optional)`:

```yaml
metadata:
  id: sampleID # str: unique identifier for the sample
  species: Bellis perennis # str: binomial name of species being assembled
  mitochondrial_code: 1 # int: mitochondrial genetic code (https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi)
sequencing_data:
  long_reads: # mandatory - information about location and platform of long reads
    platform: "pacbio" # str: "PacBio" or "ONT"
    reads: # List of str: Paths to long-read read files in FASTA or FASTQ format
      - /path/to/reads1.fa.gz
      - /path/to/reads2.fa.gz
    fastk: # (optional): if specified, all sub-fields must be present.
      directory: /path/to/fk/dir # str: path to directory containing FastK .ktab and .hist files for the long reads
      basename: basename # str: Basename of the FastK database in the directory
      kmer_size: 31 # int: Kmer size used to generate the FastK database
    coverage: 25 # (optional): int: haploid coverage of the genome to assemble.
  hic: # (optional): if not specified, no scaffolding will be run.
    reads: # List of str: Paths to Hi-C read files in CRAM format
      - /path/to/reads1.cram
      - /path/to/reads2.cram
    hic_motif: # List of str: Cutting motifs of enzymes used in Hi-C library preparation
      - GATC
      - GANTC
      - CTNAG
      - TTAA
    fastk: # (optional): if specified, all sub-fields must be present.
      directory: /path/to/fk/dir # str: path to directory containing FastK .ktab and .hist files for the Hi-C reads
      basename: basename # str: Basename of the FastK database in the directory
      kmer_size: 31 # int: Kmer size used to generate the FastK database
  illumina_10x: # (optional): if not specified, no polishing can be run.
    reads: # List of str: Paths to Illumina 10X read files in Fastq format
      - /path/to/fastq1.fq.gz
      - /path/to/fastq2.fq.gz
    fastk: # (optional): if specified, all sub-fields must be present.
      directory: /path/to/fk/dir # str: path to directory containing FastK .ktab and .hist files for the long reads
      basename: basename # str: Basename of the FastK database in the directory
      kmer_size: 31 # int: Kmer size used to generate the FastK database
  trio: # (optional) - if specified, both maternal and paternal must be specified
    maternal:
      reads: # List of str: Paths to Illumina WGS read files in Fastq format
        - /path/to/fastq1.fq.gz
        - /path/to/fastq2.fq.gz
      fastk: # (optional): if specified, all sub-fields must be present.
        directory: /path/to/fk/dir # str: path to directory containing FastK .ktab and .hist files for the long reads
        basename: basename # str: Basename of the FastK database in the directory
        kmer_size: 31 # int: Kmer size used to generate the FastK database
    paternal:
      reads: # List of str: Paths to Illumina WGS read files in Fastq format
        - /path/to/fastq1.fq.gz
        - /path/to/fastq2.fq.gz
      fastk: # (optional): if specified, all sub-fields must be present.
        directory: /path/to/fk/dir # str: path to directory containing FastK .ktab and .hist files for the long reads
        basename: basename # str: Basename of the FastK database in the directory
        kmer_size: 31 # int: Kmer size used to generate the FastK database
databases:
  busco:
    lineage: bacteria_odb10 # str: BUSCO lineage name to use for QC
  oatk:
    mito_hmm: /path/to/.hmm # str: Path to oatk HMM database
    plastid_hmm: /path/to/.hmm # str: Path to oatk HMM database
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
nextflow run sanger-tol/genomeassembly --input assets/dataset.yaml --outdir <OUTDIR> -profile docker
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
input: './samplesheet.csv'
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
