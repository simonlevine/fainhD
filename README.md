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

Then, clone this repository and run as follows:

## To run
After cloning, `cd` into this repository. Change the `input.yaml` to refer to the proper dataset (as described below). Then, invoke (with, for example, two cores):
```bash
snakemake --use-conda -j2 
```

To use an HPC, an example job script is supplied at `run.job`:
```bash
sbatch run.job
```

## Using custom data

### Local

Only paired-end read are supported. If those data are available locally (e.g., if they were downloaded before-hand), put them in the `data/raw/reads` folder and specify the library name in the `input.yaml` file. For example, if the read files were "foo_R1.fastq" and "foo_R2.fastq.fastq" and specify the sample name as:

```yaml
sample: library
```

### Sequence Read Archive

We facilitate using any paired-end RNAseq dataset from SRA. Simply specify the accession number as the sample:

```yaml
sample: SRR8250770 
```

This would run the pipeline on the HPV-positive experiment under the accession number `SRR8250770`
