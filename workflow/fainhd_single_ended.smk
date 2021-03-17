import glob

configfile: "input.yaml"

rule all:
    input:
        "reports/{}.json".format(config["sample"])

def input_file_discovery(sample):
    """given a sample name specified at the top-level, i.e., at the
    command line, determine the corresponding read files. If the
    file is not found in `data/raw/reads`, assume that this is a
    SRA accession number and download"""
    fq, = glob.glob("data/raw/reads/{}.*".format(sample))
    return fq

include: "rules/download_human_genome.smk"

rule star_single_ended:
    input:
        rules.download_genome.output["completion_flag"],
        reference_genome_dir=rules.download_genome.output[0],
        fq=lambda wildcards: input_file_discovery(wildcards.sample)
    output:
        "data/interim/{sample}_nonhost.fastq",
    conda:
        "envs/star.yaml"
    threads:
        12
    script:
        "scripts/star_se.py"

rule contig_assembly_single_ended:
    input: 
        "data/interim/{sample}_nonhost.fastq",
    output:
        "data/interim/{sample}_nonhost_contigs.fasta"
    params:
        workdir="data/interim/{sample}_spades_workdir"
    conda:
        "envs/spades.yaml"
    script:
        "scripts/spades_se.py"

# all tasks downstream of contig assembly are shared
include : "rules/shared.smk"