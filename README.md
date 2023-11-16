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

1. Parse and transform the input into the required data structures.
2. Run hifiasm in original mode.
3. Produce numerical stats, BUSCO score and QV, completeness metrics, and kmer spectra for [2].
4. If <code>hifiasm_hic_on</code> option is set, run hifiasm in HiC mode.
5. If <code>hifiasm_hic_on</code> option is set, produce numerical stats, BUSCO score and QV, completeness metrics, and kmer spectra for [4].
6. Purge primary contigs from [2], count produced primary contigs as the primary assembly.
7. Take haplotigs from [6], merge with haplotigs from [2] and purge, count produced primary contigs as assembly haplotigs.
8. Produce numerical stats, BUSCO score and QV, completeness metrics, and kmer spectra for primary and haplotigs from [6] and [7].
9. If <code>polishing_on</code> option is set, map provided 10X Illumina reads to the joined primary and alternative contigs.
10. If <code>polishing_on</code> option is set, polish initial assembly based on aligment produced in [6], then separate polished primary and hapltoigs.
11. If <code>polishing_on</code> option is set, produce numerical stats, BUSCO score and QV, completeness metrics, and kmer spectra for [10].
12. Map HiC data onto primary contigs.
13. Run scaffolding for primary contigs based on results of [12].
14. Produce numerical stats, BUSCO score and QV, completeness metrics, and kmer spectra for [13].

<img src="https://github.com/sanger-tol/genomeassembly/blob/00b25ef50f548fa6f9f3a349a45343f3fd025257/genome_assembly_pipeline.png" height="1200">

## Workflow output summary

```bash
test
├── hifiasm
│   ├── baUndUnlc1.asm.a_ctg.assembly_summary
│   ├── baUndUnlc1.asm.a_ctg.fa
│   ├── baUndUnlc1.asm.a_ctg.gfa
│   ├── baUndUnlc1.asm.ec.bin
│   ├── baUndUnlc1.asm.ovlp.reverse.bin
│   ├── baUndUnlc1.asm.ovlp.source.bin
│   ├── baUndUnlc1.asm.p_ctg.assembly_summary
│   ├── baUndUnlc1.asm.p_ctg.fa
│   ├── baUndUnlc1.asm.p_ctg.gfa
│   ├── baUndUnlc1.asm.p_utg.gfa
│   ├── baUndUnlc1.asm.r_utg.gfa
│   ├── baUndUnlc1.p_ctg.bacteria_odb10.busco
│   │   ├── baUndUnlc1-bacteria_odb10-busco
│   │   ├── baUndUnlc1-bacteria_odb10-busco.batch_summary.txt
│   │   ├── short_summary.specific.bacteria_odb10.baUndUnlc1.asm.p_ctg.fa.json
│   │   └── short_summary.specific.bacteria_odb10.baUndUnlc1.asm.p_ctg.fa.txt
│   ├── baUndUnlc1.p_ctg.ccs.merquryk
│   │   ├── baUndUnlc1.baUndUnlc1.asm.a_ctg_only.bed
│   │   ├── baUndUnlc1.baUndUnlc1.asm.a_ctg.qv
│   │   ├── baUndUnlc1.baUndUnlc1.asm.a_ctg.spectra-cn.fl.png
│   │   ├── baUndUnlc1.baUndUnlc1.asm.a_ctg.spectra-cn.ln.png
│   │   ├── baUndUnlc1.baUndUnlc1.asm.a_ctg.spectra-cn.st.png
│   │   ├── baUndUnlc1.baUndUnlc1.asm.p_ctg_only.bed
│   │   ├── baUndUnlc1.baUndUnlc1.asm.p_ctg.qv
│   │   ├── baUndUnlc1.baUndUnlc1.asm.p_ctg.spectra-cn.fl.png
│   │   ├── baUndUnlc1.baUndUnlc1.asm.p_ctg.spectra-cn.ln.png
│   │   ├── baUndUnlc1.baUndUnlc1.asm.p_ctg.spectra-cn.st.png
│   │   ├── baUndUnlc1.completeness.stats
│   │   ├── baUndUnlc1.qv
│   │   ├── baUndUnlc1.spectra-asm.fl.png
│   │   ├── baUndUnlc1.spectra-asm.ln.png
│   │   └── baUndUnlc1.spectra-asm.st.png
│   ├── polishing
│   │   ├── baUndUnlc1
│   │   ├── baUndUnlc1.consensus.fa
│   │   ├── baUndUnlc1.polished.bacteria_odb10.busco
│   │   ├── baUndUnlc1.polished.ccs.merquryk
│   │   ├── chunks
│   │   ├── haplotigs.assembly_summary
│   │   ├── haplotigs.fa
│   │   ├── merged.vcf.gz
│   │   ├── merged.vcf.gz.tbi
│   │   ├── primary.assembly_summary
│   │   ├── primary.fa
│   │   ├── refdata-baUndUnlc1
│   │   └── vcf
│   ├── purging
│   │   ├── baUndUnlc1.purged.bacteria_odb10.busco
│   │   ├── baUndUnlc1.purged.ccs.merquryk
│   │   ├── coverage
│   │   ├── coverage.htigs
│   │   ├── purged.assembly_summary
│   │   ├── purged.fa
│   │   ├── purged.htigs.assembly_summary
│   │   ├── purged.htigs.fa
│   │   ├── purge_dups
│   │   ├── purge_dups.htigs
│   │   ├── seqs
│   │   ├── seqs.htigs
│   │   ├── split_aln
│   │   └── split_aln.htigs
│   └── scaffolding
│       ├── baUndUnlc1.baUndUnlc1_scaffolds_final.ccs.merquryk
│       ├── baUndUnlc1.cram.crai
│       ├── baUndUnlc1.flagstat
│       ├── baUndUnlc1.idxstats
│       ├── baUndUnlc1.markdup.bam
│       ├── baUndUnlc1.sorted.bed
│       ├── baUndUnlc1.stats
│       └── yahs
├── hifiasm-hic
│   ├── baUndUnlc1.asm.ec.bin
│   ├── baUndUnlc1.asm.hic.a_ctg.assembly_summary
│   ├── baUndUnlc1.asm.hic.a_ctg.fa
│   ├── baUndUnlc1.asm.hic.a_ctg.gfa
│   ├── baUndUnlc1.asm.hic.hap1.p_ctg.gfa
│   ├── baUndUnlc1.asm.hic.hap2.p_ctg.gfa
│   ├── baUndUnlc1.asm.hic.p_ctg.assembly_summary
│   ├── baUndUnlc1.asm.hic.p_ctg.fa
│   ├── baUndUnlc1.asm.hic.p_ctg.gfa
│   ├── baUndUnlc1.asm.hic.p_utg.gfa
│   ├── baUndUnlc1.asm.hic.r_utg.gfa
│   ├── baUndUnlc1.asm.ovlp.reverse.bin
│   └── baUndUnlc1.asm.ovlp.source.bin
├── kmer
│   ├── baUndUnlc1_fk.hist
│   ├── baUndUnlc1_fk.ktab
│   ├── baUndUnlc1.hist
│   ├── baUndUnlc1_linear_plot.png
│   ├── baUndUnlc1_log_plot.png
│   ├── baUndUnlc1_model.txt
│   ├── baUndUnlc1_summary.txt
│   ├── baUndUnlc1_transformed_linear_plot.png
│   └── baUndUnlc1_transformed_log_plot.png
└── pipeline_info
    ├── execution_report_2023-05-24_16-09-38.html
    ├── execution_timeline_2023-05-24_16-09-38.html
    ├── execution_trace_2023-05-24_16-00-56.txt
    ├── execution_trace_2023-05-24_16-09-38.txt
    ├── pipeline_dag_2023-05-24_16-09-38.html
    └── software_versions.yml
```

## Subworkflows input summary

PREPARE_INPUT</p>

- <code>ch_input</code> - [YAML file](input_yaml) with definition of the dataset: <code>channel: datafile(yaml)</code>

GENOMESCOPE_MODEL</p>

- <code>reads</code> - Paths to reads: <code>channel: [ val(meta), [ datafile(path) ] ]</code>

RAW_ASSEMBLY</p>

- <code>hifi_reads</code> - List of files containing paths to HiFi reads: <code>channel: [ val(meta), [ datafile(path) ] ]</code>
- <code>hic_reads</code> - List of files containing paths to HiC reads: <code>channel: [ datafile(cram) ]</code>
- <code>hifiasm_hic_on</code> - Switch on/off HiC mode <code>val: Boolean</code>

PURGE_DUPS

- <code>reads_plus_assembly_ch</code> - Paths to HiFi reads, primary asm, haplotigs, genomescope model: <code>channel: [ val(meta), [ datafile(reads) ], [ datafile(pri), datafile(alt) ], datafile(model) ]</code>
- <code>prefix</code> - prefix for the output files: <code>channel: val(prefix)</code>

POLISHING

- <code>fasta_in</code> - Assembly in FASTA format with index file: <code>channel: [ val(meta), datafile(fasta), datafile(fai) ]</code>
- <code>reads_10X</code> - Path to folder with Illumina 10X FASTQ files and indices: <code>channel: datafile(path) </code>
- <code>bed_chunks_polishing</code> - Number of chunks to split fasta: <code>val: Int</code>

ALIGN_SHORT

- <code>fasta</code> - Primary contigs: <code>channel: [ val(meta), datafile(fasta) ]</code>
- <code>reads</code> - HiC reads in CRAM format <code>channel: [ val(meta), [ datafile(cram) ] ]</code>

SCAFFOLDING

- <code>bed_in</code> - Alignments coordinates after markdup: <code>channel: [ val(meta), datafile(bed) ]</code>
- <code>fasta_in</code> - Alignments coordinates after markdup: <code>channel: datafile(fasta) </code>
- <code>cool_bin</code> - Bin size for cooler: <code>val(cool_bin)</code>

GENOME_STATISTICS

- <code>assembly</code> - Primary contigs and haplotigs(optional): <code>channel: [ val(meta), datafile(pri), datafile(alt) ]</code>
- <code>lineage</code> - Paths to BUSCO database (optional) and name of the BUSCO dataset: <code>channel: [ val(meta), datafile(path), val(lineage) ]</code>
- <code>hist</code> - FASTK .hist file: <code>channel: [ val(meta), datafile(hist) ]</code>
- <code>ktab</code> - FASTK .ktab file: <code>channel: [ val(meta), datafile(ktab) ]</code>

## Subworkflows output summary

PREPARE_INPUT</p>

- <code>hifi</code> - paths to HiFi reads
- <code>hic</code> - paths to HiC reads, meta contains read group used for reaad mapping
- <code>illumina_10X</code> - paths to the folder with 10X reads
- <code>busco</code> - path to the BUSCO database (optional), name of the ODB lineage (bacteria_odb10)
- <code>primary_asm</code> - primary assembly and its indices in case provided (currently not used down the pipeline)
- <code>haplotigs_asm</code> - haplotigs assembly and its indices in case provided (currently not used down the pipeline)

GENOMESCOPE_MODEL</p>

- <code>model</code> - genomescope model
- <code>hist</code> - FASTK kmer histogram
- <code>ktab</code> - FASTK kmer table

RAW_ASSEMBLY</p>

- <code>raw_unitigs</code> - hifiasm raw unitigs in GFA format
- <code>source_overlaps</code> - hifiasm binary database of overlaps
- <code>reverse_overlaps</code> - hifiasm binary database of RC overlaps
- <code>corrected_reads</code> - hifiasm binary database of corrected reads
- <code>primary_contigs_gfa</code> - hifiasm primary contigs in GFA format
- <code>alternate_contigs_gfa</code> - hifiasm haplotigs in GFA format
- <code>processed_unitigs</code> - hifiasm processed unitigs in GFA format
- <code>primary_hic_contigs_gfa</code> - hifiasm contigs produced with integration of HiC data (in GFA format)
- <code>alternate_hic_contigs_gfa</code> - hifiasm haplotigs produced with integration of HiC data (in GFA format)
- <code>phased_hic_contigs_hap1_gfa</code> - fully phased first haplotype
- <code>phased_hic_contigs_hap2_gfa</code> - fully phased another haplotype
- <code>primary_contigs</code> - primary contigs <code>primary_contigs_gfa</code> in FASTA format
- <code>alternate_contigs</code> - haplotigs <code>alternate_contigs_gfa</code> in FASTA format
- <code>primary_hic_contigs</code> - primary contigs <code>primary_hic_contigs_gfa</code> in FASTA format
- <code>alternate_hic_contigs</code> - haplotigs <code>alternate_hic_contigs_gfa</code> in FASTA format

[Hifiasm documentation](https://hifiasm.readthedocs.io/en/latest/interpreting-output.html) contains more details

PURGE_DUPS

- <code>pri</code> - purged primary contigs
- <code>alt</code> - purged haplotigs

POLISHING

- <code>fasta</code> - polished contigs
- <code>versions</code> - versions for software used in analysis

ALIGN_SHORT

- <code>bed</code> bed file of alignments after merging and markduplicates
- <code>cram</code> cram representation of alignments
- <code>crai</code> index for <code>cram</code>
- <code>stats</code> output of samtools stats
- <code>idxstats</code> output of samtools stats idxstats
- <code>flagstat</code> output of samtools stats flagstat

SCAFFOLDING

- <code>alignments_sorted</code> output of JUICER_PRE - text file of pairs of alignmnents coordinates suitable for JUICER
- <code>fasta</code> - final scaffolds
- <code>chrom_sizes</code> - sizes of scaffolds
- <code>cool</code> - path to .cool file from COOLER_CLOAD
- <code>mcool</code> - path to .mcool file from COOLER_ZOOMIFY
- <code>snapshots</code> - image of pretext map
- <code>hic</code> - contact map in .hic format
- <code>versions</code> - versions for software used in analysis

GENOME_STATISTICS

- <code>busco</code> - busco summary in json format
- <code>merquryk_completeness</code> text file of merqury completeness score
- <code>merquryk_qv</code> text file of merqury qv score
- <code>assembly_stats_pri</code> assembly stats for the primary assembly
- <code>assembly_stats_alt</code> assembly stats for haplotigs (if provided)
- <code>versions</code> - versions for software used in analysis

## Quick Start

1. Install [`Nextflow`](https://www.nextflow.io/docs/latest/getstarted.html#installation) (`>=21.10.3`)

2. Install any of [`Docker`](https://docs.docker.com/engine/installation/), [`Singularity`](https://www.sylabs.io/guides/3.0/user-guide/) (you can follow [this tutorial](https://singularity-tutorial.github.io/01-installation/)), [`Podman`](https://podman.io/), [`Shifter`](https://nersc.gitlab.io/development/shifter/how-to-use/) or [`Charliecloud`](https://hpc.github.io/charliecloud/) for full pipeline reproducibility _(you can use [`Conda`](https://conda.io/miniconda.html) both to install Nextflow itself and also to manage software within pipelines. Please only use it within pipelines as a last resort; see [docs](https://nf-co.re/usage/configuration#basic-configuration-profiles))_.

3. Download the pipeline and test it on a minimal dataset with a single command:

   ```console
   nextflow run sanger-tol/genomeassembly -profile test,YOURPROFILE --outdir <OUTDIR>
   ```

   Note that some form of configuration will be needed so that Nextflow knows how to fetch the required software. This is usually done in the form of a config profile (`YOURPROFILE` in the example command above). You can chain multiple config profiles in a comma-separated string.

   > - The pipeline comes with config profiles called `docker`, `singularity`, `podman`, `shifter`, `charliecloud` and `conda` which instruct the pipeline to use the named tool for software management. For example, `-profile test,docker`.
   > - Please check [nf-core/configs](https://github.com/nf-core/configs#documentation) to see if a custom config file to run nf-core pipelines already exists for your Institute. If so, you can simply use `-profile <institute>` in your command. This will enable either `docker` or `singularity` and set the appropriate execution settings for your local compute environment.
   > - If you are using `singularity`, please use the [`nf-core download`](https://nf-co.re/tools/#downloading-pipelines-for-offline-use) command to download images first, before running the pipeline. Setting the [`NXF_SINGULARITY_CACHEDIR` or `singularity.cacheDir`](https://www.nextflow.io/docs/latest/singularity.html?#singularity-docker-hub) Nextflow options enables you to store and re-use the images from a central location for future pipeline runs.
   > - If you are using `conda`, it is highly recommended to use the [`NXF_CONDA_CACHEDIR` or `conda.cacheDir`](https://www.nextflow.io/docs/latest/conda.html) settings to store the environments in a central location for future pipeline runs.

4. Start running your own analysis!

   <!-- TODO nf-core: Update the example "typical command" below used to run the pipeline -->

   ```console
   nextflow run sanger-tol/genomeassembly --input sample_dataset.yaml --outdir <OUTDIR> --polishing_on true -profile <docker/singularity/podman/shifter/charliecloud/conda/institute>
   ```

## Credits

sanger-tol/genomeassembly was originally written by @ksenia-krasheninnikova.

We thank the following people for their extensive assistance in the development of this pipeline:

@priyanka-surana for nextflow implementation of HiC mapping pipeline, extensive guidance, code review, and brilliant suggestions.

@mcshane and @c-zhou for designing and implementing original pipelines for purging (@mcshane), polishing (@mcshane) and scaffolding (@c-zhou)

@mahesh-panchal for nextflow implementation of the purging pipeline, code review and the valuable suggestions for the input subworkflow

@muffato for code review and suggestions about versioning

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
