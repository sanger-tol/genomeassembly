# sanger-tol/genomeassembly: Output

## Introduction

This document describes the output produced by the pipeline.

The directories listed below will be created in the results directory after the pipeline has finished. All paths are relative to the top-level results directory.

## Kmer pre-processing subworkflow

If no pre-computed [FastK](https://github.com/thegenemyers/FASTK) databases are provided, the pipeline will generate them. These databases
are output in the `kmer` subdirectory. FastK databases will be generated for the input long reads, and if provided, the maternal and paternal
reads. For maternal and paternal reads, FastK databases are also generated containing only kmers from each parent using MerquryFK's
[hapmaker](https://github.com/thegenemyers/MERQURY.FK) tool.

If no coverage value is provided for the long reads, then [GenomeScope 2.0](https://github.com/tbenavi1/genomescope2.0) is run to estimate
the coverage of the genome prior to assembly.

<details markdown="1">
<summary>Output files</summary>

- `kmer/k{kmer_size}/long/`
  - `{id}.long.k{kmer_size}_fk.hist`: FastK Histogram file containing a histogram of kmer counts at each coverage value.
  - `{id}.long.k{kmer_size}.hist`: ASCII TSV of the FastK TSV, binned from (1..1000).
  - `{id}.long.k{kmer_size}_fk.ktab`: FastK master ktab file describing the counts of each kmer.
  - `.{id}.long.k{kmer_size}_fk.ktab.{n}`: Hidden FastK ktab files (1..n) describing the counts of each kmer.
  - `.{id}.long.k{kmer_size}_model.txt`: GenomeScope2.0 model fit file containing the model parameters.
  - `.{id}.long.k{kmer_size}_summary.txt`: Summary file for the GenomeScope2.0 model.
  - `.{id}.long.k{kmer_size}_*.png`: GenomeScope2.0 graphical plots showing the model fit against the kmer spectrum.
- `kmer/k{kmer_size}/{mat/pat}/`
  - `{id}.{mat/pat}.k{kmer_size}_fk.hist`: FastK Histogram file containing a histogram of kmer counts at each coverage value for the maternal or paternal reads.
  - `{id}.{mat/pat}.k{kmer_size}.hist`: ASCII TSV of the FastK TSV, binned from (1..1000) for the maternal or paternal reads.
  - `{id}.{mat/pat}.k{kmer_size}_fk.ktab`: FastK master ktab file describing the counts of each kmer for the maternal or paternal reads.
  - `.{id}.{mat/pat}.k{kmer_size}_fk.ktab.{n}`: Hidden FastK ktab files (1..n) describing the counts of each kmer for the maternal or paternal reads.
  - `{id}.{mat/pat}.yak`: Yak kmer database for the maternal or paternal reads.
- `kmer/{kmer_size}/trio`
  - `{id}.mat.k{kmer_size}_fk.ktab`: FastK master ktab file describing the counts of each kmer in the child present in both the child and maternal genome.
  - `.{id}.mat.k{kmer_size}_fk.ktab.{n}`: Hidden FastK ktab files (1..n) describing the counts of each kmer in the child present in both the child and maternal genome.
  - `{id}.pat.k{kmer_size}_fk.ktab`: FastK master ktab file describing the counts of each kmer in the child present in both the child and paternal genome.
  - `.{id}.pat.k{kmer_size}_fk.ktab.{n}`: Hidden FastK ktab files (1..n) describing the counts of each kmer in the child present in both the child and paternal genome.

</details>

## Assemblies

Each assembly generated is output in a separate directory in the output directory. Each assembly directory is named using the following schema:
`{id}.hifiasm(-hic/-trio)?.{date}`, where `id` is the provided sample ID, `-hic` or `-trio` are present if the assembly is Hi-C phased or
trio-binned, respectively, and `date` is the date the Nextflow workflow was started.

The following sections all describe the contents of these assembly directories.

### Hifiasm raw assembly

The raw outputs of each [hifiasm](https://github.com/chhylp123/hifiasm) assembly of the long reads are located in the base of each assembly directory.
For full details describing hifiasm outputs, please refer to the hifiasm documentation: [https://hifiasm.readthedocs.io/en/latest/interpreting-output.html](https://hifiasm.readthedocs.io/en/latest/interpreting-output.html)

<details markdown="1">
  <summary>Output files</summary>

- `asm.p_ctg.(g)?fa.gz`: primary assembly in GFA and FASTA format.
- `asm.a_ctg.(g)?fa.gz`: haplotig assembly in GFA and FASTA format.
- `asm.hap1.p_ctg.(g)fa.gz`: fully phased hap1 assembly if hifiasm is run in Hi-C phasing mode; partially phased hap1 assembly otherwise, in GFA or FASTA format.
- `asm.hap2.p_ctg.(g)fa.gz`: fully phased hap1 assembly if hifiasm is run in Hi-C phasing mode; partially phased hap1 assembly otherwise, in GFA or FASTA format.
- `asm.hap2.p_utg.gfa.gz`: haplotype-resolved processed unitig graph without small bubbles in GFA format.
- `asm.hap2.r_utg.gfa.gz`: haplotype-resolved raw unitig graph in GFA format.
- `asm.stderr.log`: hifiasm run log file.
- `*.bin`: internal binary hifiasm files. Can be used to re-run hifiasm.

</details>

### Purging

After raw assembly with hifiasm, an assembly can optionally be purged of retained haplotype using [purge_dups](https://github.com/dfguan/purge_dups). If purge_dups is run,
the purged assemblies and other associated output will be available in the `purging` directory within an assembly directory.

<details markdown="1">
  <summary>Output files</summary>

- `purging/asm.purged.fa`: purged primary contigs in FASTA format
- `purging/asm.purged.htigs.fa`: haplotigs purged from the primary assembly in FASTA format
- `purging/asm.purged.htigs.all.fa`: raw hifiasm haplotigs and haplotigs purged from the primary assembly in FASTA format
- `purging/split_aln/asm.self_aln.split.fasta.gz`: Gzipped fasta containing fragmented primary assembly for self-alignment
- `purging/coverage/asm.calcuts.log`: log file for purge_dups calcuts
- `purging/coverage/asm.cutoffs`: purge_dups cutoffs file
- `purging/coverage/asm.PB.base.cov`: purge_dups base-level read depth
- `purging/coverage/asm.PB.stat`: purge_dups read depth histogram
- `purging/coverage/asm.self_aln.paf`: PAF format self-alignment of the split primary assembly
- `purging/purge_dups/baUndUnlc1.dups.bed`: BED file describing identified retained haplotype in the primary assembly
- `purging/purge_dups/baUndUnlc1.purge_dups.log`: purge_dups log file

</details>

### Polishing

Prior to scaffolding, an assembly can be polishied using Illumina 10X reads. If purging is run, the purged assembly is polished, otherwise the raw hifiasm assembly
is polished. [Longranger](https://support.10xgenomics.com/genome-exome/software/pipelines/latest/what-is-long-ranger) and [Freebayes](https://github.com/freebayes/freebayes)
are used to polish the assemblies.

<details markdown="1">
  <summary>Output files</summary>

- `polishing/asm.consensus.fa`: Consensus polished assembly (both haplotypes) in FASTA format
- `polishing/asm_hap1.fa`: Consensus polished assembly (primary/hap1) in FASTA format
- `polishing/asm_hap2.fa`: Consensus polished assembly (alt/hap2) in FASTA format
- `polishing/asm_merged.vcf.gz`: VCF of assembly from Freebayes
- `polishing/asm_merged.vcf.gz.tbi`: TBI index of VCF of assembly from Freebayes
- `polishing/chunks/*.bed`: BED files describing assembly regions polished independently
- `polishing/asm/outs/possorted_bam.bam`: BAM file of Illumina 10X reads mapped to the combined assembly by Longranger
- `polishing/asm/outs/possorted_bam.bam.bai`: BAM index
- `polishing/asm/outs/summary.csv`: Longranger summary information

</details>

### Hi-C mapping and Hi-C mapping statistics

Illumina Hi-C reads are mapped to each assembly using either [bwa-mem2](https://github.com/bwa-mem2/bwa-mem2) or [minimap2](https://github.com/lh3/minimap2/). The mapping is performed in chunks to parallelise the process,
and the chunked BAM files are sorted by co-ordinate. The chunked BAM files are then merged, and duplicates are marked with [samtools](https://www.htslib.org/) markdup. Mapping
statistics are calculated for each mapped BAM using samtools flagstat, idxstats, and stats.

<details markdown="1">
  <summary>Output files</summary>

- `scaffolding_{hap}/asm_{hap}.bam`: Coordinate-sorted BAM file of Hi-C reads mapped to the assembly, with duplicates marked.
- `scaffolding_{hap}/asm_{hap}.flagstat`: Samtools flagstats for the BAM file of Hi-C reads mapped to the assembly.
- `scaffolding_{hap}/asm_{hap}.idxstats`: Samtools idxstats for the BAM file of Hi-C reads mapped to the assembly.
- `scaffolding_{hap}/asm_{hap}.stats`: Samtools stats for the BAM file of Hi-C reads mapped to the assembly.

</details>

### Scaffolding

Scaffolding is performed using the long-range information from the Hi-C alignments using [YaHS](https://github.com/c-zhou/yahs).
Alignments are converted to a name-sorted BED file, and then the assembly is scaffolded. Following scaffolding, Hi-C contact maps
in [Pretext](https://github.com/sanger-tol/PretextMap), [Juicer](https://github.com/aidenlab/juicer),
and [Cooler](https://github.com/open2c/cooler) are generated for visualisation. A PNG image of the Pretext map is generated for
quick visualisation of the results.

<details markdown="1">
  <summary>Output files</summary>

- `scaffolding_{hap}/asm_{hap}.bed`: Read name sorted BED file of Hi-C reads mapped to the assembly, generated by bedtools bamToBed.
- `scaffolding_{hap}/yahs/out.{break/nobreak}.yahs/asm_{hap}_scaffolds_final.fa`: final scaffolds in FASTA format
- `scaffolding_{hap}/yahs/out.{break/nobreak}.yahs/asm_{hap}_scaffolds_final.agp`: AGP file translating input contigs to scaffolds
- `scaffolding_{hap}/yahs/out.{break/nobreak}.yahs/asm_{hap}_scaffolds_final.bin`: YaHS bin file containing Hi-C contacts
- `scaffolding_{hap}/yahs/out.{break/nobreak}.yahs/asm_{hap}_scaffolds_final.pretext`: Hi-C contact map in Pretext format
- `scaffolding_{hap}/yahs/out.{break/nobreak}.yahs/asm_{hap}_scaffolds_finalFullMap.png`: PNG image of Pretext contact map
- `scaffolding_{hap}/yahs/out.{break/nobreak}.yahs/asm_{hap}_scaffolds_final.cool`: Hi-C contact map in Cooler format
- `scaffolding_{hap}/yahs/out.{break/nobreak}.yahs/asm_{hap}_scaffolds_final.hic`: Hi-C contact map in Juicer format

</details>

### Genome statistics

Accompanying every genome assembly FASTA file from the previous step are a range of genome statistics to assess the
contiguity, quality and completeness of the genome assembly. The statistics generated are: basic assembly
statistics with asmstats and [gfastats](https://github.com/vgl-hub/gfastats), BUSCO ortholog scoring with
[BUSCO](https://busco.ezlab.org/), and QV, kmer completeness and graphical visualisations
with [MerquryFK](https://github.com/thegenemyers/MERQURY.FK).

<details markdown="1">
  <summary>Output files</summary>

- `{assembly fasta}.assembly_summary`: GFAStats output for the individual assembly
- `{assembly fasta}.stats`: asmstats output for the individual assembly
- `{assembly fasta}.{busco_lineage}.busco/*`: BUSCO output directory for the individual assembly
- `{assembly fasta}.ccs.merquryk/*`: MerquryFK output directory for the pair (pri/alt, hap1/hap2) of assemblies

</details>

### Organelle assembly

Assembly of organelle genomes is accomplished with [MitoHiFi](github.com/marcelauliano/MitoHiFi) and [oatk](github.com/c-zhou/oatk).
MitoHiFi uses a reference-based approach, downloading a reference mitochondial genome for a closely-related species and
either identifying mitochondrial contigs from an existing hifiasm assembly, or assembling a mitochondial genome using reads
mapping to that reference. Oatk is a de-novo assembler that can assemble both mitochondrial and plastid genomes,
using the presence of genes to pick a path through the assembly graph.

The output of MitoHiFi in reads mode appears in the assembly directory of the standard hifiasm assembly in the directory `mito.reads`, and
the output of oatk in the directory `mito.oatk`. The output of MitoHiFi in contigs mode appears in the relevant assembly directory for the contigs
in directory `mito`.

<details markdown="1">
  <summary>Output files</summary>

- `mito(.reads)?/`: Output directory of Mitohifi in contigs mode. See the [MitoHiFi](https://github.com/marcelauliano/MitoHiFi) documentation for a full account of all files.
- `mito(.reads)?/contigs_stats.tsv`: Statistics of all mitochondrial contigs identified by MitoHiFi.
- `mito(.reads)?/final_mitogenome.fasta`: Mitochondrial genome in FASTA format chosen by MitoHiFi.
- `mito(.reads)?/final_mitogenome.gb`: Mitochondrial genome annotation in GB format chosen by MitoHiFi.
- `mito.oatk/asm.k{oatk_kmer_size}_c{oatk_coverage}.utg.final.gfa`: the GFA file for the final oatk genome assembly:
- `mito.oatk/asm.k{oatk_kmer_size}_c{oatk_coverage}.annot_mito.txt`: the MT gene annotation file for assembled sequences
- `mito.oatk/asm.k{oatk_kmer_size}_c{oatk_coverage}.annot_pltd.txt`: the PT gene annotation file for assembled sequences
- `mito.oatk/asm.k{oatk_kmer_size}_c{oatk_coverage}.mito.gfa`: the subgraph for the MT genome
- `mito.oatk/asm.k{oatk_kmer_size}_c{oatk_coverage}.mito.bed`: the gene annotation for the MT sequences
- `mito.oatk/asm.k{oatk_kmer_size}_c{oatk_coverage}.mito.ctg.fasta`: the structure-solved MT contigs
- `mito.oatk/asm.k{oatk_kmer_size}_c{oatk_coverage}.mito.ctg.bed`: the genome annotation for MT contigs
- `mito.oatk/asm.k{oatk_kmer_size}_c{oatk_coverage}.pltd.gfa`: the subgraph for the: PT genome
- `mito.oatk/asm.k{oatk_kmer_size}_c{oatk_coverage}.pltd.bed`: the gene annotation for the PT sequences
- `mito.oatk/asm.k{oatk_kmer_size}_c{oatk_coverage}.pltd.ctg.fasta`: the structure-solved PT contigs
- `mito.oatk/asm.k{oatk_kmer_size}_c{oatk_coverage}.pltd.ctg.bed`: the genome annotation for PT: contigs

</details>

## Pipeline overview

The pipeline is built using [Nextflow](https://www.nextflow.io/) and processes data using the following steps:

- [Pipeline information](#pipeline-information) - Report metrics generated during the workflow execution

### Pipeline information

<details markdown="1">
<summary>Output files</summary>

- `pipeline_info/`
  - Reports generated by Nextflow: `execution_report.html`, `execution_timeline.html`, `execution_trace.txt` and `pipeline_dag.dot`/`pipeline_dag.svg`.
  - Reports generated by the pipeline: `pipeline_report.html`, `pipeline_report.txt` and `software_versions.yml`. The `pipeline_report*` files will only be present if the `--email` / `--email_on_fail` parameter's are used when running the pipeline.
  - Reformatted samplesheet files used as input to the pipeline: `samplesheet.valid.csv`.
  - Parameters used by the pipeline run: `params.json`.

</details>

[Nextflow](https://www.nextflow.io/docs/latest/tracing.html) provides excellent functionality for generating various reports relevant to the running and execution of the pipeline. This will allow you to troubleshoot errors with the running of the pipeline, and also provide you with other information such as launch commands, run times and resource usage.
