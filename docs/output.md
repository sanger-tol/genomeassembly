# sanger-tol/genomeassembly: Output

## Introduction

This document describes the output produced by the pipeline.

The directories listed below will be created in the results directory after the pipeline has finished. All paths are relative to the top-level results directory.

Each assembly generated is output in a separate directory in the output directory. Each assembly directory is named according to the assembly specification name, specified in the `id` field of an assembly specification entry. This will be denoted `asm.id` in all documemtation that follows.

For each assembly, index JSON files are also created that describe the data and parameter inputs for each stage, as well as the software versions for all tools used. These are available in the root directory for each assembly, and additionally for nuclear assemblies there is one describing each stage individually.

## FastK Databases

If no pre-computed [FastK](https://github.com/thegenemyers/FASTK) databases for long read data typesare provided,
the pipeline will generate them. These databases are organized by sample ID, platform, and kmer size.

<details markdown="1">
<summary>Output files</summary>

- `fastk/{id}/{platform}/k{kmer_size}/`
  - `*.hist`: FastK Histogram file containing a histogram of kmer counts at each coverage value.
  - `*.ktab`: FastK master ktab file describing the counts of each kmer.
  - `.*.ktab.{n}`: Hidden FastK ktab files (1..n) describing the counts of each kmer.

</details>

> **Note:** FastK database output is only enabled when `--save_fastk_databases` is set to `true`.

## Nuclear Genome Assemblies

The following sections describe the contents of nuclear assembly directories.

### Hifiasm Raw Assembly

The raw outputs of each [hifiasm](https://github.com/chhylp123/hifiasm) assembly of the long reads are located in the `raw/` subdirectory of each assembly directory.
For full details describing hifiasm outputs, please refer to the hifiasm documentation: [https://hifiasm.readthedocs.io/en/latest/interpreting-output.html](https://hifiasm.readthedocs.io/en/latest/interpreting-output.html)

<details markdown="1">
<summary>Output files</summary>

- `{asm.id}/hifiasm/`
  - `*.fa.gz`: Assembly FASTA files (bgzip compressed). Includes primary (`p_ctg`), alternate (`a_ctg`), and/or haplotype-resolved (`hap1.p_ctg`, `hap2.p_ctg`) assemblies depending on the assembly mode.
  - `*.gfa.gz`: Assembly graphs in GFA format (bgzip compressed). Includes raw unitig graphs (`r_utg`), processed unitig graphs (`p_utg`), and primary contig graphs.
  - `*.bed.gz`: BED files with assembly information (bgzip compressed).
  - `*.log`: Hifiasm run log file.
  - `index.json`: Index file summarizing the data and parameter inputs for the stage, as well as the software versions for all tool used.

</details>

### Purging

After raw assembly with hifiasm, an assembly can optionally be purged of retained haplotype using [purge_dups](https://github.com/dfguan/purge_dups). If purge_dups is run,
the purged assemblies and other associated output will be available in the `purging/` directory within an assembly directory.

<details markdown="1">
<summary>Output files</summary>

- `{asm.id}/purging/`
  - `*.fa.gz`: Purged assembly FASTA files (bgzip compressed), including purged primary contigs and haplotigs.
- `{asm.id}/purging/split_aln/`
  - `*.split.fasta.gz`: Gzipped FASTA containing fragmented primary assembly for self-alignment.
  - `*.paf`: PAF format self-alignment of the split primary assembly.
- `{asm.id}/purging/coverage/`
  - `*.stat`: purge_dups read depth histogram (PB.stat).
  - `*.cov`: purge_dups base-level read depth (PB.base.cov).
  - `*.cutoffs`: purge_dups cutoffs file from calcuts.
  - `*.log`: Log file for purge_dups calcuts.
  - `*.png`: Histogram plot of coverage distribution.
- `{asm.id}/purging/purge_dups/`
  - `*.bed`: BED file describing identified retained haplotype in the primary assembly.
  - `*.log`: purge_dups log file.
  - `index.json`: Index file summarizing the data and parameter inputs for the stage, as well as the software versions for all tool used.

</details>

### Polishing

Prior to scaffolding, an assembly can be polished using Illumina 10X reads. If purging is run, the purged assembly is polished, otherwise the raw hifiasm assembly
is polished. [Longranger](https://support.10xgenomics.com/genome-exome/software/pipelines/latest/what-is-long-ranger) and [Freebayes](https://github.com/freebayes/freebayes)
are used to polish the assemblies.

> **Note:** Polishing requires a Longranger container to be provided via `--polishing_longranger_container_path`.

<details markdown="1">
<summary>Output files</summary>

- `{asm.id}/polishing/`
  - `*.fa.gz`: Consensus polished assembly FASTA files (bgzip compressed).
- `{asm.id}/polishing/aln/`
  - `*.bam`: BAM file of Illumina 10X reads mapped to the combined assembly by Longranger.
  - `*.bam.bai`: BAM index file.
  - `*.csv`: Longranger summary information.
- `{asm.id}/polishing/vcf/`
  - `*.vcf.gz`: Merged VCF of assembly variants from Freebayes.
  - `*.vcf.gz.tbi`: TBI index of the merged VCF.
  - `index.json`: Index file summarizing the data and parameter inputs for the stage, as well as the software versions for all tool used.

</details>

### Scaffolding

Scaffolding is performed using the long-range information from the Hi-C alignments using [YaHS](https://github.com/c-zhou/yahs).
Hi-C reads are mapped to each assembly using either [bwa-mem2](https://github.com/bwa-mem2/bwa-mem2) or [minimap2](https://github.com/lh3/minimap2/).
The mapping is performed in chunks to parallelise the process, and the chunked BAM files are sorted by coordinate, merged, and duplicates are marked with [samtools](https://www.htslib.org/) markdup.
Mapping statistics are calculated using samtools flagstat, idxstats, and stats.

Following scaffolding, Hi-C contact maps in [Pretext](https://github.com/sanger-tol/PretextMap), [Juicer](https://github.com/aidenlab/juicer),
and [Cooler](https://github.com/open2c/cooler) formats are generated for visualisation.

<details markdown="1">
<summary>Output files</summary>

- `{asm.id}/scaffolding/`
  - `*_scaffolds_final.fa.gz`: Final scaffolds in FASTA format (bgzip compressed).
- `{asm.id}/scaffolding/hic_aln/`
  - `*.bam`: Coordinate-sorted BAM file of Hi-C reads mapped to the assembly, with duplicates marked.
  - `*.csi`: CSI index for the BAM file.
  - `*.stats`: Samtools stats for the BAM file.
  - `*.flagstat`: Samtools flagstats for the BAM file.
  - `*.idxstats`: Samtools idxstats for the BAM file.
- `{asm.id}/scaffolding/yahs/`
  - `*.agp`: AGP file translating input contigs to scaffolds.
  - `*.bin`: YaHS bin file containing Hi-C contacts.
  - `*_inital.*`: Initial YaHS scaffolding outputs (if available).
  - `*_intermediate.*`: Intermediate YaHS scaffolding outputs (if available).
  - `*.log`: YaHS log file.
- `{asm.id}/scaffolding/contact_maps/`
  - `*.pretext`: Hi-C contact map in Pretext format.
  - `*FullMap.png`: PNG image of the Pretext contact map for quick visualisation.
  - `*.cool`: Hi-C contact map in Cooler format.
  - `*.hic`: Hi-C contact map in Juicer format.
  - `index.json`: Index file summarizing the data and parameter inputs for the stage, as well as the software versions for all tool used.

</details>

### Genome Statistics

Accompanying every genome assembly FASTA file are a range of genome metrics to assess the
contiguity, quality and completeness of the genome assembly. The statistics generated are: basic assembly
statistics with [asmstats](https://github.com/sanger-tol/asmstats) and [gfastats](https://github.com/vgl-hub/gfastats), ortholog completeness scoring with [BUSCO](https://busco.ezlab.org/), and QV, kmer completeness and graphical visualisations
with [MerquryFK](https://github.com/thegenemyers/MERQURY.FK).

Statistics are generated at each assembly stage (raw, purging, polishing, scaffolding) and organized in subdirectories accordingly.

<details markdown="1">
<summary>Output files</summary>

- `{asm.id}/{assembly_stage}/`
  - `*.stats`: asmstats output for the individual assembly.
  - `*.assembly_summary`: gfastats output for the individual assembly.
- `{asm.id}/{assembly_stage}/merqury.{platform}/`
  - `*.qv`: MerquryFK QV (quality value) scores.
  - `*.completeness.stats`: MerquryFK kmer completeness statistics.
  - `*.phased_block.stats`: MerquryFK phased block statistics (for trio assemblies).
  - `*.spectra-cn.fl.png`, `*.spectra-cn.ln.png`, `*.spectra-cn.st.png`: Copy number spectrum plots.
  - `*.spectra-asm.fl.png`, `*.spectra-asm.ln.png`, `*.spectra-asm.st.png`: Assembly spectrum plots.
- `{asm.id}/{assembly_stage}/busco.{lineage}/`
  - `batch_summary.txt`: BUSCO batch summary for both haplotypes.
  - `short_summary.*.txt`: BUSCO short summary text files.
  - `short_summary.*.json`: BUSCO short summary JSON files.
  - `logs/`: BUSCO log files.

</details>

### MitoHiFi

MitoHiFi uses a reference-based approach, downloading a reference mitochondrial (or plastid) genome for a closely-related species and
identifying organellar contigs from an existing hifiasm assembly. When run in contigs mode as part of a nuclear assembly pipeline,
outputs appear in the `mito/` or `plastid/` subdirectory.

<details markdown="1">
<summary>Output files</summary>

- `{asm.id}/{organelle}/`
  - `final_mitogenome.fasta`: Final mitochondrial/plastid genome in FASTA format.
  - `final_mitogenome.gb`: Mitochondrial/plastid genome annotation in GenBank format.
  - `final_mitogenome.gff`: Mitochondrial/plastid genome annotation in GFF format.
  - `contigs_stats.tsv`: Statistics of all organellar contigs identified.
  - `all_potential_contigs.fa`: All potential organellar contigs identified.
  - `contigs_annotations.png`: Visualization of contig annotations.
  - `coverage_plot.png`: Coverage plot across the organellar genome.
  - `final_mitogenome.annotation.png`: Annotation visualization of the final genome.
  - `final_mitogenome.coverage.png`: Coverage visualization of the final genome.
  - `shared_genes.tsv`: Shared genes information.
  - `*.log`: MitoHiFi log file.
  - `contigs_circularization/`: Directory containing circularization information.
  - `contigs_filtering/`: Directory containing filtering information.
  - `coverage_mapping/`: Directory containing coverage mapping data.
  - `final_mitogenome_choice/`: Directory containing information about final genome selection.
- `{asm.id}/{organelle}/ref/`
  - `*.fa`: Reference organellar genome FASTA used by MitoHiFi.
  - `*.gb`: Reference organellar genome GenBank annotation.
  - `index.json`: Index file summarizing the data and parameter inputs for the stage, as well as the software versions for all tool used.

</details>

## Organelle Assembly

Assembly of organelle genomes is accomplished with [MitoHiFi](https://github.com/marcelauliano/MitoHiFi) and [oatk](https://github.com/c-zhou/oatk).

### MitoHiFi

When MitoHiFi is run in reads mode, reads are matched to the downloaded reference assembly and assembled. The
outputs are largely the same as when MitoHifi is run on assembled contigs.

<details markdown="1">
<summary>Output files</summary>

- `{asm.id}/`
  - `final_mitogenome.fasta`: Final mitochondrial/plastid genome in FASTA format.
  - `final_mitogenome.gb`: Mitochondrial/plastid genome annotation in GenBank format.
  - `final_mitogenome.gff`: Mitochondrial/plastid genome annotation in GFF format.
  - `contigs_stats.tsv`: Statistics of all organellar contigs identified.
  - `all_potential_contigs.fa`: All potential organellar contigs identified.
  - `contigs_annotations.png`: Visualization of contig annotations.
  - `coverage_plot.png`: Coverage plot across the organellar genome.
  - `final_mitogenome.annotation.png`: Annotation visualization of the final genome.
  - `final_mitogenome.coverage.png`: Coverage visualization of the final genome.
  - `shared_genes.tsv`: Shared genes information.
  - `*.log`: MitoHiFi log file.
  - `contigs_circularization/`: Directory containing circularization information.
  - `contigs_filtering/`: Directory containing filtering information.
  - `coverage_mapping/`: Directory containing coverage mapping data.
  - `final_mitogenome_choice/`: Directory containing information about final genome selection.
- `{asm.id}/ref/`
  - `*.fa`: Reference organellar genome FASTA used by MitoHiFi.
  - `*.gb`: Reference organellar genome GenBank annotation.
  - `index.json`: Index file summarizing the data and parameter inputs for the stage, as well as the software versions for all tool used.

</details>

### Oatk

[Oatk](https://github.com/c-zhou/oatk) is a de-novo assembler that can assemble both mitochondrial and plastid genomes,
using the presence of genes to pick a path through the assembly graph.

<details markdown="1">
<summary>Output files</summary>

- `{asm.id}/`
  - `*.mito.ctg.fasta`: Structure-solved mitochondrial contigs in FASTA format.
  - `*.pltd.ctg.fasta`: Structure-solved plastid contigs in FASTA format.
  - `*.mito.ctg.bed`: Genome annotation BED file for mitochondrial contigs.
  - `*.pltd.ctg.bed`: Genome annotation BED file for plastid contigs.
  - `*.mito.gfa`: Subgraph for the mitochondrial genome in GFA format.
  - `*.pltd.gfa`: Subgraph for the plastid genome in GFA format.
  - `*.annot_mito.txt`: Mitochondrial gene annotation file for assembled sequences.
  - `*.annot_pltd.txt`: Plastid gene annotation file for assembled sequences.
  - `*.utg.gfa`: Initial oatk unitig graph.
  - `*.utg.final.gfa`: Final oatk unitig graph.
  - `*.log`: Oatk log file.
  - `*.png`: Bandage plots visualizing the assembly graphs.
  - `index.json`: Index file summarizing the data and parameter inputs for the stage, as well as the software versions for all tool used.

</details>

## Pipeline Information

The pipeline is built using [Nextflow](https://www.nextflow.io/) and processes data using the following steps:

<details markdown="1">
<summary>Output files</summary>

- `pipeline_info/`
  - `execution_report.html`: Nextflow execution report.
  - `execution_timeline.html`: Nextflow execution timeline.
  - `execution_trace.txt`: Nextflow execution trace.
  - `pipeline_dag.dot`, `pipeline_dag.svg`: Pipeline directed acyclic graph visualizations.
  - `genomeassembly_software_versions.yml`: Software versions used in the pipeline run.
  - `pipeline_report.html`, `pipeline_report.txt`: Pipeline reports (only present if `--email` / `--email_on_fail` parameters are used).
  - `params.json`: Parameters used by the pipeline run.

</details>

[Nextflow](https://docs.seqera.io/platform-cloud/reports/overview) provides excellent functionality for generating various reports relevant to the running and execution of the pipeline. This will allow you to troubleshoot errors with the running of the pipeline, and also provide you with other information such as launch commands, run times and resource usage.
