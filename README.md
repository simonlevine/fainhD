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

   - STAR Alignment
    - User supplies the reference genome sequences (FASTA files) and annotations (GTF file), from which STAR generate genome indexes.
    - By default, we use and **currently only support a precomputed STAR index file**.
    - The genome indexes need only be generated once for each reference genome/annotation combination.
  
  - Mapping reads: In this step user supplies the genome files generated in the 1st step, as well as the RNA-seq reads (sequences) in the form of FASTA or FASTQ files. STAR maps the reads to the genome, and writes several output files, such as alignments (SAM/BAM), mapping summary statistics, splice junctions, unmapped reads, signal (wiggle) tracks etc.

  - Extract Unmapped Reads
    - We aim to isolate non-host sequences from the RNA-seq data.
    - We will use this to use contig

2. Assemble unknown sequences into contigs
   - In progress...
3. BLAST virus sequences against known viruses
   - In progress...
4. Predict functional ORFs in viral sequences
   - In progress...
5. Search for structural elements in virus sequences
   - In progress...