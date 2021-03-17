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

rule star_se:
    input:
        rules.download_genome.output["completion_flag"],
        fq1=lambda wildcards: input_file_discovery(wildcards.sample)
    output:
        unmapped_reads="data/interim/star/{sample}/Unmapped.out.mate1"
    params:
        index="data/raw/reference_genome",
        extra="--outReadsUnmapped Fastx"
    conda:
        "envs/star.yaml"
    wrapper:
        "0.72.0/bio/star/align"

rule contig_assembly_se:
    input: 
        rules.star_se.output["unmapped_reads"]
    output:
        "data/interim/{sample}_nonhost_contigs.fasta"
    params:
        sequencing="single-ended"
    threads:
        12
    conda:
        "envs/spades.yaml"
    script:
        "scripts/spades.py"

# all tasks downstream of contig assembly are shared
include : "rules/shared.smk"