[![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.XXXXXXX-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.XXXXXXX)

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A521.10.3-23aa62.svg)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![Launch on Nextflow Tower](https://img.shields.io/badge/Launch%20%F0%9F%9A%80-Nextflow%20Tower-%234256e7)](https://tower.nf/launch?pipeline=https://github.com/sanger-tol/genomeassembly)

## Introduction

<!-- TODO nf-core: Write a 1-2 sentence summary of what data the pipeline is for and what it does -->

**sanger-tol/genomeassembly** is a bioinformatics pipeline for a genome assembly for HiFi, Illumina 10x (optional), and HiC data. It performs the following steps: raw assembly, purging from haplotigs, optional polishing, and scaffolding.

Original assembly of HiFi reads is performed using [hifiasm](https://hifiasm.readthedocs.io) assembler in two modes - original and using HiC data (optional). Then assembly is purged from alternative haplotigs using [purge_dups](https://github.com/dfguan/purge_dups). Next optional step is polishing of the purged assembly using Illumina 10X read sequencing. 10X reads are mapped to the full assembly (purged + haplotigs) using [Longranger](https://support.10xgenomics.com/genome-exome/software/pipelines/latest/what-is-long-ranger) and polishing is implemented using [Freebayes](https://github.com/freebayes/freebayes). HiC reads are further mapped with [bwamem2](https://github.com/bwa-mem2/bwa-mem2) to the primary contigs, which are further scaffolded with [YaHS](https://github.com/c-zhou/yahs) using the provided Hi-C data.
Polished and scaffolded assemblies are evaluated using [GFASTATS](https://github.com/vgl-hub/gfastats), [BUSCO](https://busco.ezlab.org/) and [MERQURY.FK](https://github.com/thegenemyers/MERQURY.FK)

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It uses Docker/Singularity containers making installation trivial and results highly reproducible. The [Nextflow DSL2](https://www.nextflow.io/docs/latest/dsl2.html) implementation of this pipeline uses one container per process which makes it much easier to maintain and update software dependencies. Where possible, these processes have been submitted to and installed from [nf-core/modules](https://github.com/nf-core/modules) in order to make them available to all nf-core pipelines, and to everyone within the Nextflow community!

<!-- TODO nf-core: Add full-sized test dataset and amend the paragraph below if applicable -->

On release, automated continuous integration tests run the pipeline on a full-sized dataset on the AWS cloud infrastructure. This ensures that the pipeline runs on AWS, has sensible resource allocation defaults set to run on real-world datasets, and permits the persistent storage of results to benchmark between pipeline releases and other analysis sources.

## Pipeline summary

While the steps are described in a sequential order, many of them can be executed as parallel jobs.

1. Parse and transform the input into the required data structures.
2. Run organnels subworkflow on the HiFi reads.
3. Run hifiasm in the original mode.
4. Produce numerical stats, BUSCO score and QV, completeness metrics, and kmer spectra for [3].
5. If <code>hifiasm_hic_on</code> option is set, run hifiasm in HiC mode.
6. If <code>hifiasm_hic_on</code> option is set, produce numerical stats, BUSCO score and QV, completeness metrics, and kmer spectra for [5].
7. Purge primary contigs from [3], i.e. produce the purged assembly and a set of haplotigs. Consider the purged contigs as the primary assembly for further steps.
8. Take haplotigs from [7], merge with haplotigs from [3] and purge again. Discard the purhed contigs, continue with the purged haplotigs as representation of the alternative assembly further on.
9. If not <code>polishing_on</code> run organelles subwrokflow on the joined purged contigs and haplotigs.
10. Produce numerical stats, BUSCO score and QV, completeness metrics, and kmer spectra for primary and haplotigs from [8] and [9].
11. If <code>polishing_on</code> option is set, map Illumina 10X reads to the joined primary and alternative contigs.
12. If <code>polishing_on</code> option is set, polish initial assembly based on aligment produced in [8], then separate polished primary and hapltoigs.
13. If <code>polishing_on</code> option is set, produce numerical stats, BUSCO score and QV, completeness metrics, and kmer spectra for [12].
14. If <code>polishing_on</code> run organelles subwrokflow on the joined polished contigs and haplotigs.
15. Map HiC data onto primary contigs.
16. Run scaffolding for primary contigs based on results of [15].
17. Produce numerical stats, BUSCO score and QV, completeness metrics, and kmer spectra for [16].

## Quick Start

1. Install [`Nextflow`](https://www.nextflow.io/docs/latest/getstarted.html#installation) (`>=21.10.3`)

2. Install any of [`Docker`](https://docs.docker.com/engine/installation/), [`Singularity`](https://www.sylabs.io/guides/3.0/user-guide/) (you can follow [this tutorial](https://singularity-tutorial.github.io/01-installation/)), [`Podman`](https://podman.io/), [`Shifter`](https://nersc.gitlab.io/development/shifter/how-to-use/) or [`Charliecloud`](https://hpc.github.io/charliecloud/) for full pipeline reproducibility _(you can use [`Conda`](https://conda.io/miniconda.html) both to install Nextflow itself and also to manage software within pipelines. Please only use it within pipelines as a last resort; see [docs](https://nf-co.re/usage/configuration#basic-configuration-profiles))_.

3. Download the pipeline and test it on a minimal dataset with a single command:

   ```console
   nextflow run sanger-tol/genomeassembly -profile test,YOURPROFILE --outdir <OUTDIR>
   ```

   Note that some form of configuration will be needed so that Nextflow knows how to fetch the required software. This is usually done in the form of a config profile (`YOURPROFILE` in the example command above). You can chain multiple config profiles in a comma-separated string.

Further documentation about the pipeline can be found in the following files: usage, parameters and output.

4. Start running your own analysis!

   <!-- TODO nf-core: Update the example "typical command" below used to run the pipeline -->

   ```console
   nextflow run sanger-tol/genomeassembly --input sample_dataset.yaml --outdir <OUTDIR> --polishing_on true -profile <docker/singularity/podman/shifter/charliecloud/conda/institute>
   ```

## Credits

sanger-tol/genomeassembly was originally written by @ksenia-krasheninnikova.

We thank the following people for their extensive assistance in the development of this pipeline:

@priyanka-surana for the code review, brilliant coding suggestions, and help with pushing this pipeline forward through development.

@mcshane and @c-zhou for designing and implementing original pipelines for purging (@mcshane), polishing (@mcshane) and scaffolding (@c-zhou).

@muffato for help with nf-core code and troubleshooting, code review and valuable suggestions at the different stages of pipeline development.

@mahesh-panchal for nextflow implementation of the purging pipeline, code review and valuable suggestions to the nf-core modules implementation.

<!-- TODO nf-core: If applicable, make list of people who have also contributed -->

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

## Citations

<!-- TODO nf-core: Add citation for pipeline after first release. Uncomment lines below and update Zenodo doi and badge at the top of this file. -->
<!-- If you use  sanger-tol/genomeassembly for your analysis, please cite it using the following doi: [10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX) -->

<!-- TODO nf-core: Add bibliography of tools and data used in your pipeline -->

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

This pipeline uses code and infrastructure developed and maintained by the [nf-core](https://nf-co.re) community, reused here under the [MIT license](https://github.com/nf-core/tools/blob/master/LICENSE).

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
