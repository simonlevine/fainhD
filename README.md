![](fainhD.png)

**f**iltering **a**nd **i**dentifying **n**on-**h**ost **D**NA pipeline for 03-713 Bioinformatics Practicum at Carnegie Mellon University, Spring 2021 by Jeremy Fisher, Simon Levine, Sid Reed and Tomas Matteson

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

- [x] Filter RNA-seq data to remove host sequences
- [x] Assemble unknown sequences into contigs
- [ ] BLAST virus sequences against known viruses
- [ ] Predict functional ORFs in viral sequences
- [ ] Search for structural elements in virus sequences
