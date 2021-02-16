
# rule star_pe_multi: #paired-end sequence
#     input:
#         # use a list for multiple fastq files for one sample
#         # usually technical replicates across lanes/flowcells
#         fq1 = ["./input/alignment/reads/{sample}_R1.1.fastq", "./input/alignmentreads/{sample}_R1.2.fastq"],
#         # # paired end reads needs to be ordered so each item in the two lists match
#         fq2 = ["./input/alignment/reads/{sample}_R2.1.fastq", "./input/alignmentreads/{sample}_R2.2.fastq"] #optional
#     output:
#         # see STAR manual for additional output files
#         "star/pe/sample_sequence/Aligned.out.sam"
#     log:
#         "logs/star/pe/{sample}.log"
#     params:
#         # path to STAR reference genome index
#         index="./input/alignment/index",
#         # optional parameters
#         extra="--outSAMunmapped Within"
#         # extra="--foo bar"

#     threads: 8
#     wrapper:
#         "0.72.0/bio/star/align"

rule star_se: #single-end
    input:
        # fq1 = "./input/alignment/reads/{sample}_R1.1.fastq"
        fq1 = 'data/sample_sequence.fq'

    output:
        # see STAR manual for additional output files
        "./output/alignment/sample/Aligned.out.sam" #alignments in standard SAM format

    log:
        "./output/logs/star/sample.log"
    params:
        # path to STAR reference genome index
        index="index",
        # optional parameters
        # extra=""
        extra="--outSAMunmapped Within"

    threads: 8
    wrapper:
        "0.72.0/bio/star/align"



# rule samtools_sort:
#     input:
#         "mapped_reads/{sample}.bam"
#     output:
#         "sorted_reads/{sample}.bam"
#     shell:
#         "samtools sort -T sorted_reads/{wildcards.sample} "
#         "-O bam {input} > {output}"


# rule samtools_index:
#     input:
#         "sorted_reads/{sample}.bam"
#     output:
#         "sorted_reads/{sample}.bam.bai"
#     shell:
#         "samtools index {input}"


# rule bcftools_call:
#     input:
#         fa="data/genome.fa",
#         bam=expand("sorted_reads/{sample}.bam", sample=SAMPLES),
#         bai=expand("sorted_reads/{sample}.bam.bai", sample=SAMPLES)
#     output:
#         "calls/all.vcf"
#     shell:
#         "samtools mpileup -g -f {input.fa} {input.bam} | "
#         "bcftools call -mv - > {output}"


# rule plot_quals:
#     input:
#         "calls/all.vcf"
#     output:
#         "plots/quals.svg"
#     script:
#         "scripts/plot-quals.py"