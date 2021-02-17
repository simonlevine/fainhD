from snakemake.io import directory
from snakemake.remote.HTTP import RemoteProvider

HTTP = RemoteProvider()

reference_genome_url_prefix = "http://labshare.cshl.edu/shares/gingeraslab/www-data/dobin/STAR/STARgenomes/Human/GRCh38_Ensembl99_sparseD3_sjdbOverhang99"

rule all:
    input: "data/interim/sample_aligned2humangenome.sam"

rule download_genome:
    input:
        [HTTP.remote(f"{reference_genome_url_prefix}/{f}", keep_local=True) for f in ['chrLength.txt', 'chrName.txt', 'chrStart.txt', 'Genome', 'genomeParameters.txt', 'SA', 'SAindex']]
    output:
        outdir=directory("data/raw/reference_genome"),
        download_complete_flag="data/raw/reference_genome/DOWNLOAD_COMPLETE.txt"
    run:
        for f in input:
            shell("mv {f} {output.outdir}")
        shell("touch {output.download_complete_flag}")

rule star_single_ended:
    input: 
        "data/raw/reference_genome/DOWNLOAD_COMPLETE.txt",
        fq1 = "data/raw/{sample}_sequence.fq"
    output:
        "data/interim/{sample}_aligned2humangenome.sam"
    log:
        "logs/star/{sample}.log"
    params:
        # path to STAR reference genome index
        index="data/raw/reference_genome",
        extra="--outSAMunmapped Within"
    threads: 8
    wrapper:
        "0.72.0/bio/star/align"

#rule extract_unmapped_reads:
#    input:
#        "./output/alignment/{sample}/Aligned.out.sam"
#    output:
#        "./output/filtering/{sample}/nonhost_sequences.sam"
#    params:
#        # '4' is the flag for unmapped reads
#        # see: http://www.htslib.org/doc/samtools.html
#        "-f 4"
#    wrapper:
#        "0.72.0/bio/samtools/view" 
#
#rule convert_sam_to_fastq:
#    input:
#        "./output/filtering/{sample}/nonhost_sequences.sam"
#    output:
#        "./output/filtering/{sample}/nonhost_sequences.fastq"
#    shell:
#        # see: https://www.cureffi.org/2013/07/04/how-to-convert-sam-to-fastq-with-unix-command-line-tools/
#        """grep -v ^@ | awk '{print "@"$1"\n"$10"\n+\n"$11}' < {input} > {output}"""
#
#rule contig_assembly:
#    input: 
#        "./output/filtering/{sample}/nonhost_sequences.fastq"
#    output:
#        directory("./output/contigs/{sample}")
#    conda:
#        "environment.yaml"
#    shell:
#        "spades.py --rna --s1 {input} -o {output}"
