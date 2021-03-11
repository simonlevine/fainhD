# Under The Hood

## Snakemake

We employ Snakemake to build a directed acyclic graph of tasks to execute stepwise.

## RNA-Seq Filtering: STAR

### Algorithm
STAR operates on an indexed genome (specifically, a suffix array), against which subsequences of a `fasta` are queried. For a given query, the longest subsequence that matches a location in the genome is stored as a "seed". These seeds are extended, but they need not match precisely. Because mismatches are allowed, many locations in the genome may "align" to a given query. All such alignments that are allowed are scored and, ultimately, only the best alignment is recorded.

### Snakemake Rule
```
rule star_double_ended:
    input:
        rules.download_genome.output["completion_flag"],
        reference_genome_dir=rules.download_genome.output[0],
        fq1 = lambda wildcards: input_file_validation(wildcards.sample, "1"),
        fq2 = lambda wildcards: input_file_validation(wildcards.sample, "2")
    output:
        "data/interim/{sample}_nonhost_R1.fastq",
        "data/interim/{sample}_nonhost_R2.fastq"
    conda:
        "envs/star.yaml"
    threads:
        12
    script:
        "scripts/star.py"
```

## Contig Assembly: SPADES

### Algorithm

We use rnaSPAdes based on the original SPAdes:

1. Assembly graph construction (multi-sized de Bruijn graph)
   i. aggregates biread information into distance histograms, among others.
2. *k*-bimer adjustment
   i. derives accurate distance estimates between k-mers in the genome (edges in the assembly graph) using joint analysis of such distance histograms
3. Construct paired assembly graph
4. Contig construction
   i. Construct DNA sequences of contigs and the mapping of reads to contigs

Further Reading:
- https://academic.oup.com/gigascience/article/8/9/giz100/5559527
- https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3342519/#s035

### Snakemake Rule

```
rule contig_assembly:
    input: 
        "data/interim/{sample}_nonhost_R1.fastq",
        "data/interim/{sample}_nonhost_R2.fastq"
    output:
        "data/interim/{sample}_nonhost_contigs.fasta"
    params:
        workdir="data/interim/{sample}_spades_workdir"
    conda:
        "envs/spades.yaml"
    script:
        "scripts/spades.py"
```

## Alignment to Known Viruses: BLASTn
After you have isolated reads from your data you suspect are viruses you want to be able to see if your suspected viral reads are part of a known virus.
So a BLAST search of your suspected viral contigs against all known human viruses is preformed so you can see if any of your contigs are similar to known viruses.
If not it may be that a novel virus is present in your data or the contig is not actually a virus but an human or transposons sequence (among other possibilities).
It is also possible your virus has certain similar genes or motifs to a known virus but is still different enough to be considered novel.

### Algorithm
The input is a query sequence and a set of search sequences called the database, the output is the similarity of the input to every sequence in the output and associated statistics.
To grossly oversimplify, the query sequence is chopped up into equally sized substrings which can individually be matched against sequences in the database.

BLASTn uses a heuristic version of the Smith-Waterman algorithm for local sequence alignment, intended to increase the speed of searching queries against a database.
1. Filter low-complexity regions from the query sequence
1. Chop query sequence into k-mers
1. Match and score each unique k-mer againt the database
1. Only keep matched k-mers with a score above the threshold
1. Exactly match each kept k-mer against the database
1. Preform a local alignment in both directions around each exact match k-mer
1. Keep any high scoring local alignments and statistically validate those matches

Further Reading
- [Having a BLAST with bioinformatics (and avoiding BLASTphemy)](https://genomebiology.biomedcentral.com/articles/10.1186/gb-2001-2-10-reviews2002#Sec7)
- [wikipedia](https://en.wikipedia.org/wiki/BLAST)
- [video overview](https://www.youtube.com/watch?v=LlnMtI2Sg4g)
- [original paper](https://www.sciencedirect.com/science/article/pii/S0022283605803602?via%3Dihub)

### Snakemake Rule
First a local blast database of all known viruses must be created.
So a fasta file of all known __complete__ viral genomes is downloaded from NCBI (specifically [this file](https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.1.1.genomic.fna.gz)).
This file is downloaded and then the command `makeblastdb` is used to make a blast database (the 5 `*.n*` files in the `ext` list).
Now any query can be blasted against this database, which contains all known complete human viral genomes (~9000 as of March 2021).

```
rule make_known_virus_blastdb:
    input:
        "data/raw/viral_genomes.fasta"
    params:
        db_prefix="data/blast/virus_blastdb",  # used in downstream rule, should be a global?
        virus_fasta='data/raw/viral_genomes.fasta'
    output:
        expand("data/blast/virus_blastdb.{ext}", ext= ["nhr", "nin", "nog", "nsq"])
    conda:
        "envs/blast.yaml"
    shell:
        "makeblastdb -dbtype nucl "
        "-parse_seqids "
        "-input_type fasta "
        "-in {params.virus_fasta} "
        "-out {params.db_prefix} "
```

Now that the blast database is made we can blast each individual contig against the database to see which viruses it is most similar to using Smith-Waterman local alignment.
The `-outfmt 6` flag is described in the `blastn -h` help page, but it produces a tsv file with every hit from the given contig against the database.
Note that there may be many low quality hits since `blastn` is not filtering the output, this is left to the user to parse as different users will be looking for different answers in their blast results.

```
rule blast_against_known_viruses:
    input:
        "data/blast/virus_blastdb.nhr",
        contigs="data/interim/{sample}_nonhost_contigs.fasta"
    output:
        "data/processed/{sample}_blast_results.tsv"
    params:
        db_name="data/blast/virus_blastdb" # createde upstream, should be global?
    conda:
        "envs/blast.yaml"
    shell:
        "blastn -db {params.db_name} "
        "-outfmt 6 "
        "-query {input.contigs} "
        "-out {output}"
```

The .tsv output of `blastn` is parsed in to an easier to parse .csv for each input sample to `data/processed/{sample}_blast_results.tsv`.
The most important columns are
- `qseqid` The fasta header of the suspected viral contig input
- `sseqid` The fasta header of the known virus it matched to
- `pident` The percent identity of the query (qseqid) that matched the subject (sseqid)
- `evalue` Statistical evaluation of whether this match is spurious, see [here](https://www.ncbi.nlm.nih.gov/BLAST/tutorial/Altschul-1.html#head4) for details

You can also look online or on the `blastn` man page for more details about the outputs and their meaning.

```
rule convert_blast_to_csv_with_header:
    input:
        "data/processed/{sample}_blast_results.tsv"
    output:
        "data/processed/{sample}_blast_results.csv"
    run:
        cols = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
                "qstart", "qend", "sstart", "send", "evalue", "bitscore"]
        with open(input[0]) as fin, open(output[0], "wt") as fout:
            tsv_reader = csv.reader(fin, delimiter="\t")
            csv_writer = csv.DictWriter(fout, fieldnames=cols)
            csv_writer.writeheader()
            for row in tsv_reader:
                csv_writer.writerow(dict(zip(cols, row)))
```

If you want all results in a single file you can simply run the following in the main project directory.
```
$ cat data/processed/*_blast_results.csv >| all_blast_results.csv
```

## Functional ORF Prediction: Orfipy

### Algorithm
Orfipy is based on Aho–Corasick string-searching:

1. Construct a trie out of input patterns.
2. Construct suffix and output links in Breadth First Order.
3. Use the constructed automata to traverse the string.
4. If current character matches any of children then follow it otherwise follow the suffix link.
5. At every node follow the output links to get patterns occurring till the current position.

Complexity:
- Trie Contruction: $O(n)$
- Suffix/Output Link Construction: $O(n)$
- Searching: $O(m + z)$
- Time Complexity: $O(n + m + z)$

Futher Reading:
- https://github.com/urmi-21/orfipy
- https://iq.opengenus.org/aho-corasick-algorithm

### Snakemake Rule

We output a final `.fasta` file per the Snakemake rule:
```
rule pORF_finding:
    input:
        "data/interim/{sample}_nonhost_contigs.fasta"
    output:
        "data/processed/{sample}_orf.fasta"
    conda:
        "envs/orfipy.yaml"
    shell:
        "orfipy {input} --rna {output} --min 10 --max 10000 --table 1 --outdir ."

```
### Algorithm

We use rnaSPAdes based on the original SPAdes:

1. Assembly graph construction (multi-sized de Bruijn graph)
   i. aggregates biread information into distance histograms, among others.
2. *k*-bimer adjustment
   i. derives accurate distance estimates between k-mers in the genome (edges in the assembly graph) using joint analysis of such distance histograms
3. Construct paired assembly graph
4. Contig construction
   i. Construct DNA sequences of contigs and the mapping of reads to contigs

Further Reading:
- https://academic.oup.com/gigascience/article/8/9/giz100/5559527
- https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3342519/#s035

### Snakemake Rule

```
rule contig_assembly:
    input:
        "data/interim/{sample}_nonhost_R1.fastq",
        "data/interim/{sample}_nonhost_R2.fastq"
    output:
        "data/interim/{sample}_nonhost_contigs.fasta"
    params:
        workdir="data/interim/{sample}_spades_workdir"
    conda:
        "envs/spades.yaml"
    script:
        "scripts/spades.py"
```


## Structural Inference: Inferal

### Algorithm

Infernal works by querying a covariance model database and performing multiple sequence alignments between the models in that database and your input sequence. Each discovered structural RNA family has 1 structural covariance model, generated from many seed alignments that can all be found on Rfam.

Rfam has around 3500 structural RNA families, and more are constantly being added (which can also be done by building CMs from sequence on Infernal). A covariance model is a specialized type of stochastic contact free grammar, related to a profile Hidden Markov Model. The difference being (at a very high level), that in profile HMMs each nucleotide at each position is independent, but in CMs the nucleotides are dependent upon each other.
The actual explanation is a lot more complicated, and involves knowledge of Mixture Dirichlet priors, an extended discussion on which can be found in Chapter 5 of Infernal's manual.
How do we use it?

First we download the whole Rfam database as a covariance model database, and then we run Infernal on our fasta file and cm database. For each sequence in our FASTA, the input sequence is compared against all the CMs in the database, and based on pre-calibrated features from Rfam, if a significant hit is reported, the position and identity of that hit in the sequence is returned.

Infernal considers overlapping sequences, so we can have the case where from positions 0-100 we have RNA family 1 and from 80-120 we have RNA family 2, etc. 

Infernal must match on known covariance models, as it is essientially a very complicated multiple sequence alignment, so never before seen RNA structures must be characterized and added to the database before Infernal can report finding them.

### Snakemake Rule

```
rule structural_prediction:
    input:
        query_fasta="data/interim/{sample}_viral_contigs.fasta",
        rfam_database="data/raw/rfam/Rfam.cm"
    output:
        tblout="data/processed/{sample}_structural_predictions.tblout",
        cmscan="data/processed/{sample}_structural_predictions.cmscan"
    conda:
        "envs/infernal.yaml"
    threads:
        12
    shell:
        "cmscan --rfam --cut_ga --nohmmonly "
        "--cpu {threads} "
        "--tblout {output.tblout} " 
        "{input.rfam_database} "
        "{input.query_fasta} "
        "> {output.cmscan}"
```