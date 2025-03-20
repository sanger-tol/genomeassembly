# sanger-tol/genomeassembly: Usage

## Introduction

## Workflow input

### Parameters summary

<details markdown="1">
  <summary>Details</summary>
  
Workflow accepts the following parameters:
* <code>input</code> - (required) YAML file containing description of the dataset, incl. ToLID, paths to the raw data etc.
* <code>bed_chunks_polishing</code> - a number of chunks to split contigs for polishing (default 100)
* <code>cool_bin</code> - a bin size for cooler (default 1000)
* <code>organelles_on</code> - set <code>True</code> for running organelles subworkflow
* <code>polishing_on</code> - set <code>True</code> for polishing
* <code>hifiasm_hic_on</code> - set <code>True</code> to run of hifiasm in HiC mode
  <br>NB: hifiasm in the original mode is used as the main assembly even if the <code>hifiasm_hic_on</code> flag is set</br>

</details>

### Full samplesheet

The input dataset is described in YAML format which states for "Yet Another Markdown Language". It is a human readable file which contains information
about location paths for the raw data (HiFi, 10X, HiC) used for the genome assembly. It can also contain meta information such as HiC restriction motifs,
BUSCO lineage, mitochondrial code etc. For more information see [Input YAML definition](#input_yaml_definition)

### <a name="input_yaml_definition"></a> Input YAML definition

- <code>dataset.id</code>
  - is used as the sample id throughout the pipeline. ToLID should be used in ToL datasets.</p>
- <code>dataset.illumina_10X.reads</code>
  - is necessary in case polishing is applied, this field should point to the path of the folder containing 10X reads. Sample identifier in the Illumina reads should coincide with the top level ID. For the use of the Longranger software the reads should follow [the 10X FASTQ file naming convention](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/fastq-input).</p>
- <code>dataset.pacbio.reads</code>
  - contains the list (<code>-reads</code>) of the HiFi reads in FASTA (or gzipped FASTA) format in. The pipeline implementation is based on an assumption that reads have gone through adapter/barcode checks. </p>
- <code>dataset.HiC.reads</code>
  - contains the list (<code>-reads</code>) of the HiC reads in the indexed CRAM format.</p>
- <code>dataset.hic_motif</code>
  - is a comma-separated list of restriction sites. The pipeline was tested with the Arima dataset, but it's should be alright to use it with the other HiC libraries
- <code>dataset.busco.lineage</code>
  - specifies the name of the BUSCO dataset (i.e. bacteria_odb10). </p>
- <code>dataset.busco.lineage_path</code>
  - is an optional field containing the path to the folder with pre-downloaded BUSCO lineages. </p>
- <code>dataset.mito.species</code>
  - is the latin name of the species to look for the mitogenome reference in the organelles subworkflow. Normally this parameter will contain the latin name of the species whose genome is being assembled. </p>
- <code>dataset.mito.min_length</code>
  - sets the minimal length of the mito, can be 15Kb. </p>
- <code>dataset.mito.code</code>
  - is a mitochondrial code for the mitogenome annotation. See [here](https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi) for reference.

### <a name="input_yaml_example"></a> An example of the input YAML

<details markdown="1">
  <summary>Details</summary>
  
Example is based on [test.yaml](../assets/test.yaml).
```yaml
dataset:
  id: baUndUnlc1
  illumina_10X:
    reads:
      - https://tolit.cog.sanger.ac.uk/test-data/Undibacterium_unclassified/genomic_data/baUndUnlc1/10x/baUndUnlc1_S12_L002_R1_001.fastq.gz
      - https://tolit.cog.sanger.ac.uk/test-data/Undibacterium_unclassified/genomic_data/baUndUnlc1/10x/baUndUnlc1_S12_L002_R2_001.fastq.gz
      - https://tolit.cog.sanger.ac.uk/test-data/Undibacterium_unclassified/genomic_data/baUndUnlc1/10x/baUndUnlc1_S12_L002_I1_001.fastq.gz
  pacbio:
    reads:
      - reads: https://tolit.cog.sanger.ac.uk/test-data/Undibacterium_unclassified/genomic_data/baUndUnlc1/pacbio/fasta/HiFi.reads.fasta
  HiC:
    reads:
      - reads: https://tolit.cog.sanger.ac.uk/test-data/Undibacterium_unclassified/genomic_data/baUndUnlc1/hic-arima2/41741_2%237.sub.cram
hic_motif: GATC,GANTC,CTNAG,TTAA
hic_aligner: bwamem2
busco:
  lineage: bacteria_odb10
mito:
  species: Caradrina clavipalpis
  min_length: 15000
  code: 5
  fam: https://github.com/c-zhou/OatkDB/raw/main/v20230921/insecta_mito.fam
plastid:
  fam: https://github.com/c-zhou/OatkDB/raw/main/v20230921/acrogymnospermae_pltd.fam
```
</details>

## Extra installation procedures

### Longranger

Longranger is a proprietary software product from 10X Genomics.
Its terms and conditions state that we _cannot_ redistribute the copy we have in the Tree of Life department.

If you want to run the polising option, you have to install longranger yourself.
Go to <https://support.10xgenomics.com/genome-exome/software/downloads/latest>,
read their End User Software License Agreement,
and you'll be able to download the software if you accept it.

To make a Docker (or Singularity) container out of it,
use the following Dockerfile.

```Dockerfile
FROM ubuntu:22.04
LABEL org.opencontainers.image.licenses="10x Genomics End User Software License Agreement - https://support.10xgenomics.com/genome-exome/software/downloads/latest"
ARG DEST=/opt
ADD ./longranger-2.2.2.tar.gz $DEST
RUN ln -s $DEST/longranger-2.2.2/longranger /usr/local/bin/
```

Then, to use the container in the pipeline, write the following to a `longranger.config` file

```
process {
    withName: LONGRANGER_MKREF {
        container = "/path/to/longranger_container"
    }

    withName: LONGRANGER_ALIGN {
        container = "/path/to/longranger_container"
    }
}
```

And pass it to the pipeline with `-c longranger.config`.

## Usage

### Local testing

<details markdown="1">
  <summary>Details</summary>

The pipeline can be tested locally using a provided small test dataset:

```
git clone git@github.com:sanger-tol/genomeassembly.git
cd genomeassembly/
nextflow run main.nf -profile test,singularity --outdir ${OUTDIR} {OTHER ARGUMENTS}
```

These command line steps will download the pipeline and run the test.

You should now be able to run the pipeline as you see fit.

</details>

### Running the pipeline

The typical command for running the pipeline is as follows:

```console
nextflow run sanger-tol/genomeassembly --input assets/dataset.yaml --outdir <OUTDIR> -profile docker
```

This will launch the pipeline with the `docker` configuration profile, also using your institution profille if available (see [nf-core/configs](#nf-core_configs)). See below for more information about profiles.

In case organelles subworkflow is switched on you will also need to set a nextflow secret to store the API key belonging to your user.

```bash
  nextflow secrets set TOL_API_KEY '[API key]'
```

Note that the pipeline will create the following files in your working directory:

```console
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

You can also generate such `YAML`/`JSON` files via [nf-core/launch](https://nf-co.re/launch).

### Updating the pipeline

When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. When running the pipeline after this, it will always use the cached version if available - even if the pipeline has been updated since. To make sure that you're running the latest version of the pipeline, make sure that you regularly update the cached version of the pipeline:

```console
nextflow pull sanger-tol/genomeassembly
```

### Reproducibility

It is a good idea to specify the pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.

First, go to the [sanger-tol/genomeassembly releases page](https://github.com/sanger-tol/genomeassembly/releases) and find the latest version number - numeric only (eg. `1.3.1`). Then specify this when running the pipeline with `-r` (one hyphen) - eg. `-r 1.3.1`.

This version number will be logged in reports when you run the pipeline, so that you'll know what you used when you look back in the future.

To further assist in reproducibility, you can use share and reuse [parameter files](#running-the-pipeline) to repeat pipeline runs with the same settings without having to write out a command with every single parameter.

> [!TIP]
> If you wish to share such profile (such as upload as supplementary material for academic publications), make sure to NOT include cluster specific paths to files, nor institutional specific profiles.

## Core Nextflow arguments

> [!NOTE]
> These options are part of Nextflow and use a _single_ hyphen (pipeline parameters use a double-hyphen)

### `-profile`

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments.

Several generic profiles are bundled with the pipeline which instruct the pipeline to use software packaged using different methods (Docker, Singularity, Podman, Shifter, Charliecloud, Conda) - see below. When using Biocontainers, most of these software packaging methods pull Docker containers from quay.io e.g [FastQC](https://quay.io/repository/biocontainers/fastqc) except for Singularity which directly downloads Singularity images via https hosted by the [Galaxy project](https://depot.galaxyproject.org/singularity/) and Conda which downloads and installs software locally from [Bioconda](https://bioconda.github.io/).

> [!IMPORTANT]
> We highly recommend the use of Docker or Singularity containers for full pipeline reproducibility, however when this is not possible, Conda is also supported.

The pipeline also dynamically loads configurations from [https://github.com/nf-core/configs](https://github.com/nf-core/configs) when it runs, making multiple config profiles for various institutional clusters available at run time. For more information and to check if your system is supported, please see the [nf-core/configs documentation](https://github.com/nf-core/configs#documentation).

Note that multiple profiles can be loaded, for example: `-profile test,docker` - the order of arguments is important!
They are loaded in sequence, so later profiles can overwrite earlier profiles.

If `-profile` is not specified, the pipeline will run locally and expect all software to be installed and available on the `PATH`. This is _not_ recommended, since it can lead to different results on different machines dependent on the computer environment.

- `docker`
  - A generic configuration profile to be used with [Docker](https://docker.com/)
- `singularity`
  - A generic configuration profile to be used with [Singularity](https://sylabs.io/docs/)
- `podman`
  - A generic configuration profile to be used with [Podman](https://podman.io/)
- `shifter`
  - A generic configuration profile to be used with [Shifter](https://nersc.gitlab.io/development/shifter/how-to-use/)
- `charliecloud`
  - A generic configuration profile to be used with [Charliecloud](https://hpc.github.io/charliecloud/)
- `apptainer`
  - A generic configuration profile to be used with [Apptainer](https://apptainer.org/)
- `wave`
  - A generic configuration profile to enable [Wave](https://seqera.io/wave/) containers. Use together with one of the above (requires Nextflow ` 24.03.0-edge` or later).
- `conda`
  - A generic configuration profile to be used with [Conda](https://conda.io/docs/). Please only use Conda as a last resort i.e. when it's not possible to run the pipeline with Docker, Singularity, Podman, Shifter or Charliecloud.
- `test`
  - A profile with a complete configuration for automated testing
  - Includes links to test data so needs no other parameters

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

```console
NXF_OPTS='-Xms1g -Xmx4g'
```
