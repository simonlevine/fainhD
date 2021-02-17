# 03-713 Bioinformatics Practicum, Carnegie Mellon University, Spring 2021
## Jeremy Fisher, Simon Levine, Sid Reed, Tomas Matteson
## Pipeline

1. Filter RNA-seq data to remove host sequences
   i. STAR Alignment
    - User supplies the reference genome sequences (FASTA files) and annotations (GTF file), from which STAR generate genome indexes.
    - By default, we use and **currently only support a precomputed STAR index file**.
      - Precomputed STAR index (~100mb):
        - http://labshare.cshl.edu/shares/gingeraslab/www-data/dobin/STAR/STARgenomes/Human/GRCh38_Ensembl99_sparseD3_sjdbOverhang99/
        - You'll need `chrLength.txt  chrName.txt  chrStart.txt  Genome  genomeParameters.txt  SA  SAindex`
    - The genome indexes need only be generated once for each reference genome/annotation combination.
  
  ii. Mapping reads: In this step user supplies the genome files generated in the 1st step, as well as the RNA-seq reads (sequences) in the form of FASTA or FASTQ files. STAR maps the reads to the genome, and writes several output files, such as alignments (SAM/BAM), mapping summary statistics, splice junctions, unmapped reads, signal (wiggle) tracks etc.

  iii. Extract Unmapped Reads
    - We aim to isolate non-host sequences from the RNA-seq data.
    - We will use this to use contig

2. Assemble unknown sequences into contigs
   - In progress.
3. BLAST virus sequences against known viruses
4. Predict functional ORFs in viral sequences
5. Search for structural elements in virus sequences

### Appendix:

1. Download reference index and STAR genome data
  - Run `downloadSTAR.sh`
2. STAR Indexing:
 - Example genome:
     - http://labshare.cshl.edu/shares/gingeraslab/www-data/dobin/STAR/STARgenomes/Human/GRCh38_Ensembl99_sparseD3_sjdbOverhang99/Homo_sapiens.GRCh38.dna.primary_assembly.fa
 - Example annotations:
     - http://labshare.cshl.edu/shares/gingeraslab/www-data/dobin/STAR/STARgenomes/Human/GRCh38_Ensembl99_sparseD3_sjdbOverhang99/Homo_sapiens.GRCh38.99.gtf
