![](pipeline-logo.png)

**f**iltering **a**nd **i**dentifying **n**on-**h**ost **D**NA pipeline, for 03-713 Bioinformatics Practicum at Carnegie Mellon University, Spring 2021

By Jeremy Fisher, Simon Levine, Sid Reed, Tomas Matteson

## Installation

Please refer to [Snakemake's installation instructions](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html), as Snakemake is the only dependency. (It will, in turn, handle the pipelines software and data dependencies.)

## To run:

To run with the default dataset: 

```bash
snakemake --use-conda -j1 all
```

An example invocation is supplied for HPCs, for example:

```bash
sbatch run.job
```

## Description

fainhD does the following

1. Filter RNA-seq data to remove host sequences
2. Assemble unknown sequences into contigs
3. BLAST virus sequences against known viruses
4. Predict functional ORFs in viral sequences
5. Search for structural elements in virus sequences