configfile: "config/config.yaml"

from snakemake.io import directory
from snakemake.remote.HTTP import RemoteProvider

HTTP = RemoteProvider()
reference_genome_url_prefix = config["reference_genome_url_prefix"]

rule all:
    input: "data/interim/demo_aligned_to_human_genome.sam"

rule download_genome:
    input:
        [HTTP.remote(f"{reference_genome_url_prefix}/{f}", keep_local=True)
         for f in ['chrLength.txt', 'chrName.txt', 'chrStart.txt', 'Genome', 'genomeParameters.txt', 'SA', 'SAindex']]
    output:
	directory("data/raw/reference_genome")
    run:
        for f in input:
	    shell("mv {f} {output}"


rule star_double_ended:
    input:
        fq1 = "data/raw/{sample}_R1.fastq",
        fq2 = "data/raw/{sample}_R2.fastq",
	reference_genome_dir=rules.download_genome.output[0]
    output:
        "data/interim/{sample}_aligned_to_human_genome.sam",
    conda:
        "workflow/envs/star.yaml"
    shell:
	"STAR"
	"--genomeDir {input.reference_genome_dir} "
	"--outFileNamePrefix {output}/ "

rule star_single_ended:
    input: 
        fq1 = "data/raw/{sample}.fq",
        reference_genome_dir="data/raw/reference_genome",
    output:
        "data/interim/{sample}_aligned_to_human_genome.sam"
    log:
        "logs/star/{sample}.log"
    params:
        # path to STAR reference genome index
        index="data/raw/reference_genome",
        extra="--outSAMunmapped Within"
    threads: 8
    wrapper:
        "0.72.0/bio/star/align"

rule extract_unmapped_reads:
   input:
        "data/interim/{sample}_aligned_to_human_genome.sam"
   output:
        "data/interim/{sample}_nonhost.sam"
   params:
       # '4' is the flag for unmapped reads
       # see: http://www.htslib.org/doc/samtools.html
       "-f 4"
   wrapper:
       "0.72.0/bio/samtools/view" 

rule convert_sam_to_fastq:
    input:
        "data/interim/{sample}_nonhost.sam"
    output:
        "data/interim/{sample}_nonhost.fq"
    shell:
        # see: https://www.cureffi.org/2013/07/04/how-to-convert-sam-to-fastq-with-unix-command-line-tools/
        """grep -v ^@ | awk '{print "@"$1"\n"$10"\n+\n"$11}' < {input} > {output}"""

rule contig_assembly:
    input: 
        "data/interim/{sample}_nonhost.fq"
    output:
        "data/interim/{sample}_nonhost_contigs.fasta"
    params:
        workdir="data/interim/{sample}_spades_workdir"
    conda:
        "workflow/envs/spades.yaml"
    shell:
        "spades.py --rna --s1 {input} -o {params.workdir} "
        "&& mv {params.workdir}/contigs.fasta {output} "
        "&& rm -rf {params.workdir}"
