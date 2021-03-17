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

rule star_pe:
    input:
        rules.download_genome.output["completion_flag"],
        fq1 = lambda wildcards: input_file_discovery(wildcards.sample, "1"),
        fq2 = lambda wildcards: input_file_discovery(wildcards.sample, "2")
    output:
        outdir="data/interim/star/{sample}",
        unmapped_reads_R1="data/interim/star/{sample}/Unmapped.out.mate1",
        unmapped_reads_R2="data/interim/star/{sample}/Unmapped.out.mate2"
    params:
        index="data/raw/reference_genome",
        extra="--outReadsUnmapped Fastx"
    threads:
        12
    conda:
        "envs/star.yaml"
    wrapper:
        "0.72.0/bio/star/align"

rule contig_assembly_pe:
    input: 
        rules.star_pe.output["unmapped_reads_R1"],
        rules.star_pe.output["unmapped_reads_R2"]
    output:
        "data/interim/{sample}_nonhost_contigs.fasta"
    params:
        sequencing="pair-ended"
    threads:
        12
    conda:
        "envs/spades.yaml"
    script:
        "scripts/spades.py"

# all tasks downstream of contig assembly are shared
include : "rules/shared.smk"