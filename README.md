![](fainhD.png)

**f**iltering **a**nd **i**dentifying **n**on-**h**ost **D**NA pipeline for 03-713 Bioinformatics Practicum at Carnegie Mellon University, Spring 2021 by Jeremy Fisher, Simon Levine, Sid Reed and Tomas Matteson

fainhD does the following

- Filter RNA-seq data to remove host sequences
- Assemble unknown sequences into contigs
- BLAST virus sequences against known viruses
- Predict functional ORFs in viral sequences
- Search for structural elements in virus sequences

It takes an Illumina paired-read sequencing files and report these into a `json-lines` file with the following columns:

- `query_name`: arbitrary name assigned by the contig assembly software
- `contiq_seq`: contiguous DNA sequence itself
- `rfam_e_value`: expect value for a query into the rfam structure database
- `blast_results`: paired loci and expect values from BLAST, with keys `sseqid` for the loci and `evalue` for the expect value

This report, once completed, is deposited in `report/{sample}.json`, where `{sample}` is the name of library.

Useful intermediary files are available in the `data/processed` directory, including:

- `{sample}_blast_results.csv`: comma-delimited table of BLAST alignments to determine the identity of viral sequences
- `{sample}_orf.fasta`: open reading frame predictions
- `{sample}_structural_predictions.cmscan` visualized alignments against the rfam database, demonstrating RNA secondary structure (if any is found)
- `{sample}_structural_predictions.tblout` alignments against the rfam database, summarized

## Installation

We organize the pipeline using Snakemake.
[![Snakemake](https://img.shields.io/badge/snakemake-≥5.6.0-brightgreen.svg?style=flat)](https://snakemake.readthedocs.io)

Please refer to [Snakemake's installation instructions](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html), as Snakemake is the only dependency. It will, in turn, handle the pipelines software and data dependencies (besides the `fastq` files themselves).

Note that you may need to run ```pip install ftputil``` if you plan to use remote files (as is the default).

Then, clone this repository and run as follows:

## To run
After cloning, `cd` into this repository. Change the `input.yaml` to refer to the proper dataset (as described below). Then, invoke (with, for example, twelve cores):
```bash
snakemake --use-conda -j12 
```

To use an HPC, an example job script is supplied at `run.job`:
```bash
sbatch run.job
```

## Using custom data

### Local

Only paired-end read are supported. If those data are available locally (e.g., if they were downloaded before-hand), put them in the `data/raw/reads` folder and specify the library name in the `input.yaml` file. For example, if the read files were "foo_R1.fastq" and "foo_R2.fastq.fastq" and specify the sample name as:

```yaml
sample: foo
```

### Sequence Read Archive

We facilitate using any paired-end RNAseq dataset from SRA. Simply specify the accession number as the sample:

```yaml
sample: SRR8250770 
```

This would run the pipeline on the HPV-positive experiment under the accession number `SRR8250770`

### Batch Invocation

Snakemake is built to produce target files from source files. By specifying a `sample`, internally the pipeline requests the target file `report/{sample}.json`. This, in turn, produces the other dependencies in the pipeline like, for example, the contigous reads `fasta` file.

This makes it simple to analyze multiple libraries simultaneously. For example, to analyze the SRA experiments `SRR8250770` and `SRR8250765` (which are HPV-positive and HPV-negative, respectively):

```bash
snakemake --use-conda -j12 reports/{SRR8250770,SRR8250765}.json
```

It shall be noted that this is simply a bash expansion. The shell actually runs the following 

```bash
snakemake --use-conda -j12 reports/SRR8250770.json reports/SRR8250765.json
```