# 03-713 Bioinformatics Practicum, Carnegie Mellon University, Spring 2021

By Jeremy Fisher, Simon Levine, Sid Reed, Tomas Matteson

## Installation

Please refer to [Snakemake's installation instructions](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html), as Snakemake is the only dependency. (It will, in turn, handle the pipelines software and data dependencies.)

## Invocation

To run with the default dataset: 

```bash
snakemake --use-conda -j1 all
```

## Description

1. Filter RNA-seq data to remove host sequences

   - STAR Alignment: In this step user supplies the RNA-seq reads (sequences) in the form of paired FASTQ files. STAR maps the reads to the genome, and writes several output files, such as alignments (SAM/BAM), mapping summary statistics, splice junctions, unmapped reads, signal (wiggle) tracks etc.

  - Extract Unmapped Reads
    - We aim to isolate non-host sequences from the RNA-seq data.
    - We will use this to create contiguous sequences

2. Assemble unknown sequences into contigs
   - In progress...
3. BLAST virus sequences against known viruses
   - In progress...
4. Predict functional ORFs in viral sequences
   - In progress...
5. Search for structural elements in virus sequences
   - In progress...