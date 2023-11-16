[![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.XXXXXXX-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.XXXXXXX)

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A522.10.1-23aa62.svg)](https://www.nextflow.io/)
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

1. Parse the input into channels.
2. Run organnels subworkflow on the HiFi reads.
3. Run hifiasm in the original mode.
4. Produce numerical stats, BUSCO score and QV, completeness metrics, and kmer spectra for [3].
5. If <code>hifiasm_hic_on</code> option is set, run hifiasm in HiC mode.
6. If <code>hifiasm_hic_on</code> option is set, produce numerical stats, BUSCO score and QV, completeness metrics, and kmer spectra for [5].
7. Run purging subworkflow on the primary contigs from [3], i.e. produce the purged assembly and a set of haplotigs. Consider the purged contigs as the primary assembly for further steps.
8. Take haplotigs from [7], merge with haplotigs from [3] and run purging subworkfllow on it. Discard the contigs that were purged away, continue with the purged haplotigs as a representation of the haplotig assembly.
9. Produce numerical stats, BUSCO score and QV, completeness metrics, and kmer spectra for the primary and haplotigs from [7] and [8].
10. If <code>polishing_on</code> option is set, map Illumina 10X reads to the joined primary and alt contigs.
11. If <code>polishing_on</code> option is set, polish initial assembly based on the aligment produced in [11]. Set polished primary contigs as the primary assembly and polished haplotigs as the haplotig assembly.
12. If <code>polishing_on</code> option is set, produce numerical stats, BUSCO score and QV, completeness metrics, and kmer spectra for [12].
13. Run organelles subworkflow on the joined primary and haplotigs contigs.
14. Map HiC data onto primary contigs.
15. Run scaffolding for primary contigs based on results in [15].
16. Produce numerical stats, BUSCO score and QV, completeness metrics, and kmer spectra for [16].

## Usage

> **Note**
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/usage/installation) on how
> to set-up Nextflow. Make sure to [test your setup](https://nf-co.re/docs/usage/introduction#how-to-run-a-pipeline)
> with `-profile test` before running the workflow on actual data.

Currently, it is advised to run the pipeline with docker or singularity as some modules do not have a conda env associated with them.

To run the pipeline use a command-line:

```bash
nextflow run sanger-tol/genomeassembly -profile singularity,YOURPROFILE --outdir <OUTDIR>
```

For more details on how to run the pipeline and interprete the results see [usage](https://github.com/sanger-tol/genomeassembly/blob/dev/docs/usage.md) and [output](https://github.com/sanger-tol/genomeassembly/blob/dev/docs/output.md) sections of the documentation.

## Credits

sanger-tol/genomeassembly was originally written by @ksenia-krasheninnikova.

We thank the following people for their extensive assistance in the development of this pipeline:

@priyanka-surana for the code review, very helpful coding suggestions, and assistance with pushing this pipeline forward through development.

@mcshane and @c-zhou for the design and implementation of the original pipelines for purging (@mcshane), polishing (@mcshane) and scaffolding (@c-zhou).

TreeVal team Damon-Lee Pointon (@DLBPointon), Yumi Sims (@yumisims) and William Eagles (@weaglesBio) for implementation of the hic-mapping pipeline.

@muffato for help with nf-core integration, dealing with infrastructure and troubleshooting, for the code reviews and valuable suggestions at the different stages of the pipeline development.

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
