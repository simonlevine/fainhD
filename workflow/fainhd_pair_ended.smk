import re
import csv
from typing import Literal
from pathlib import Path

configfile: "input.yaml"

def input_file_validation(sample, paired_end: Literal["1", "2"]):
    """given a sample name specified at the top-level, i.e., at the
    command line, determine the corresponding read files. If the
    file is not found in `data/raw/reads`, assume that this is a
    SRA accession number and download"""
    data_dir = Path("data/raw/reads")
    pat = r"^" + sample + r"_?R" + paired_end + r"_?\.(fastq|fq)(\.gz)?"
    fqs = [fq for fq in data_dir.glob("*") if re.match(pat, fq.name)]
    try:
        fq, = fqs
    except ValueError:
        print(f"Could not find paired end file R{paired_end} for sample {sample}. "
              f"Assuming this is a SRA accession number...")
        return str(data_dir/f"{sample}_{paired_end}.fastq")
    return str(fq)

rule all:
    input:
        "reports/{}.json".format(config["sample"])

include: "rules/download_from_sequence_read_archives.smk"
include: "rules/download_human_genome.smk"

rule star_pair_ended:
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

rule contig_assembly_pair_ended:
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

include: "rules/download_viral_genomes.smk"

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

rule pORF_finding:
    input:
        "data/interim/{sample}_viral_contigs.fasta"
    output:
        "data/processed/{sample}_orf.fasta"
    conda:
        "envs/orfipy.yaml"
    shell:
        "orfipy {input} --rna {output} --min 10 --max 10000 --table 1 --outdir ."

rule filter_out_nonhits:
    input:
        "data/processed/{sample}_blast_results.tsv",
        "data/interim/{sample}_nonhost_contigs.fasta",
    output:
        "data/interim/{sample}_viral_contigs.fasta"
    conda:
        "envs/biopython.yaml"
    script:
        "scripts/filter_out_non_blast_hits.py"

include: "rules/download_rfam.smk"

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

rule report:
    input:
        "data/processed/{sample}_structural_predictions.tblout",
        "data/processed/{sample}_blast_results.csv",
        "data/interim/{sample}_viral_contigs.fasta",
        "data/processed/{sample}_orf.fasta"
    output:
        "reports/{sample}.json"
    conda:
        "envs/biopython.yaml"
    script:
        "scripts/report.py"
    
