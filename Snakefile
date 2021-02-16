

rule star_pe_multi:
    input:
        # use a list for multiple fastq files for one sample
        # usually technical replicates across lanes/flowcells
        fq1 = ["./input/alignment/reads/{sample}_R1.1.fastq", "./input/alignmentreads/{sample}_R1.2.fastq"],
        # paired end reads needs to be ordered so each item in the two lists match
        fq2 = ["./input/alignment/reads/{sample}_R2.1.fastq", "./input/alignmentreads/{sample}_R2.2.fastq"] #optional
    output:
        # see STAR manual for additional output files
        "star/pe/{sample}/Aligned.out.sam"
    log:
        "logs/star/pe/{sample}.log"
    params:
        # path to STAR reference genome index
        index="./input/alignment/index",
        # optional parameters
        extra=""
    threads: 8
    wrapper:
        "0.72.0/bio/star/align"

rule star_se:
    input:
        fq1 = "./input/alignment/reads/{sample}_R1.1.fastq"
    output:
        # see STAR manual for additional output files
        "./output/alignment/{sample}/Aligned.out.sam" #alignments in standard SAM format
        
    log:
        "./output/logs/star/{sample}.log"
    params:
        # path to STAR reference genome index
        index="index",
        # optional parameters
        extra=""
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
