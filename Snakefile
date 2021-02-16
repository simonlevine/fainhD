rule star_single_ended:
    input:
        fq1 = 'input/alignment/reads/novel/sample_sequence.fq'

    output:
        # see STAR manual for additional output files
        "./output/alignment/sample/Aligned.out.sam" #alignments in standard SAM format
    log:
        "./output/logs/star/sample.log"
    params:
        # path to STAR reference genome index
        index="input/reference",
        extra="--outSAMunmapped Within"
    threads: 8
    wrapper:
        "0.72.0/bio/star/align"

rule extract_unmapped_reads:
    input:
        "./output/alignment/{sample}/Aligned.out.sam"
    output:
        "./output/filtering/{sample}/nonhost_sequences.sam"
    params:
        # '4' is the flag for unmapped reads
        # see: http://www.htslib.org/doc/samtools.html
        "-f 4"
    wrapper:
        "0.72.0/bio/samtools/view" 

rule convert_sam_to_fastq:
    input:
        "./output/filtering/{sample}/nonhost_sequences.sam"
    output:
        "./output/filtering/{sample}/nonhost_sequences.fastq"
    shell:
        # see: https://www.cureffi.org/2013/07/04/how-to-convert-sam-to-fastq-with-unix-command-line-tools/
        """grep -v ^@ | awk '{print "@"$1"\n"$10"\n+\n"$11}' < {input} > {output}"""

rule contig_assembly:
    input: 
        "./output/filtering/{sample}/nonhost_sequences.fastq"
    output:
        directory("./output/contigs/{sample}")
    conda:
        "environment.yaml"
    shell:
        "spades.py --rna --s1 {input} -o {output}"