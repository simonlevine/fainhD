import re
import csv
from typing import Literal
from pathlib import Path

configfile: "input.yaml"

rule all:
    input:
        "reports/{}.json".format(config["sample"])

def input_file_discovery(sample, paired_end: Literal["1", "2"]):
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

include: "rules/download_from_sequence_read_archives.smk"
include: "rules/download_human_genome.smk"

rule star_pair_ended:
    input:
        rules.download_genome.output["completion_flag"],
        reference_genome_dir=rules.download_genome.output[0],
        fq1 = lambda wildcards: input_file_discovery(wildcards.sample, "1"),
        fq2 = lambda wildcards: input_file_discovery(wildcards.sample, "2")
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

# all tasks downstream of contig assembly are shared
include : "rules/shared.smk"