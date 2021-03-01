import os
from snakemake.remote.HTTP import RemoteProvider

HTTP = RemoteProvider()
reference_genome_url_prefix = "http://labshare.cshl.edu/shares/gingeraslab/www-data/dobin/STAR/STARgenomes/Human/GRCh38_Ensembl99_sparseD3_sjdbOverhang99"

rule download_from_SRA:
    output:
        "../data/raw/reads/{accession}_1.fastq",
        "../data/raw/reads/{accession}_2.fastq",
    wrapper:
        "0.72.0/bio/sra-tools/fasterq-dump" 

rule download_genome:
    input:
        [HTTP.remote(f"{reference_genome_url_prefix}/{f}", keep_local=True)
         for f in ['chrLength.txt', 'chrName.txt', 'chrStart.txt', 'Genome',
                   'genomeParameters.txt', 'SA', 'SAindex', "sjdbInfo.txt"]]
    output:
        directory("../data/raw/reference_genome"),
        completion_flag=touch("../data/raw/reference_genome/download_genome.done")
    run:
        for f in input:
            shell("mv {f} {output[0]}")

rule star_double_ended:
    input:
        rules.download_genome.output["completion_flag"],
        fq1 = [os.path.join("..", fp) for fp in config["sequence_files_R1"]],
        fq2 = [os.path.join("..", fp) for fp in config["sequence_files_R2"]],
        reference_genome_dir=rules.download_genome.output[0]
    output:
        "../data/interim/{sample}_nonhost_R1.fastq",
        "../data/interim/{sample}_nonhost_R2.fastq"
    conda:
        "../envs/star.yaml"
    threads: 12
    shell:
        "export STAR_WORK_DIR=../data/interim/{wildcards.sample}_alignment_workingdir/ "
        "&& mkdir -p $STAR_WORK_DIR "
        "&& STAR"
        " --runThreadN {threads}"
        " --genomeDir {input.reference_genome_dir}"
        " --readFilesIn {input.fq1} {input.fq2}"
        " --outFileNamePrefix $STAR_WORK_DIR"
        " --outReadsUnmapped Fastx "
        "&& mv $STAR_WORK_DIR/Unmapped.out.mate1 {output[0]} "
        "&& mv $STAR_WORK_DIR/Unmapped.out.mate2 {output[1]} "
        ";  rm -rf $STAR_WORK_DIR"