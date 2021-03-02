![](fainhD.png)

**f**iltering **a**nd **i**dentifying **n**on-**h**ost **D**NA pipeline for 03-713 Bioinformatics Practicum at Carnegie Mellon University, Spring 2021 by Jeremy Fisher, Simon Levine, Sid Reed and Tomas Matteson

## Description

fainhD does the following

- [x] Filter RNA-seq data to remove host sequences
- [x] Assemble unknown sequences into contigs
- [ ] BLAST virus sequences against known viruses
- [x] Predict functional ORFs in viral sequences
- [ ] Search for structural elements in virus sequences

## Installation

Please refer to [Snakemake's installation instructions](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html), as Snakemake is the only dependency. It will, in turn, handle the pipelines software and data dependencies.

Then, clone this repository and 

## To run

After cloning, `cd` into this repository. Change the `input.yaml` to refer to the proper dataset (as described below). Then, invoke:
```bash
snakemake --use-conda -j1 all
```

To use an HPC, an example job script is supplied at `run.job`:
```bash
sbatch run.job
```

## Using custom data

### Local

Only paired-end read data are supported. If those data are available locally (e.g., if they were downloaded before-hand), specify their filepaths in `input.yaml` relative to the repository root:

```yaml
sequence_files_R1: 
- "data/raw/reads/library_R1.fastq"
sequence_files_R2: 
- "data/raw/reads/library_R2.fastq"
```

In this case, we expect read libraries to reside in `data/raw/reads`, although they can reside anywhere 

### Sequence Read Archive

We facilitate using any paired-end RNAseq dataset from SRA. Simply specify the accession number like so in the `input.yaml`:
```yaml
sequence_files_R1: 
- "data/raw/reads/{accession-number}_1.fastq"
sequence_files_R2: 
- "data/raw/reads/{accession-number}_2.fastq"
```
For example, to use the HPV-positive experiment under `SRR8250770`, use:
```yaml
sequence_files_R1: 
- "data/raw/reads/SRR8250770_1.fastq"
sequence_files_R2: 
- "data/raw/reads/SRR8250770_2.fastq"
```

Please note, if one does use the automatic SRA download feature, the `data/raw/reads` bit is **not optional**.
