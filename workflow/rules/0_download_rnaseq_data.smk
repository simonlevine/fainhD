from snakemake.remote.HTTP import RemoteProvider
HTTP = RemoteProvider()

rule download_from_SRA:
    output:
        "../data/raw/reads/{accession}_1.fastq",
        "../data/raw/reads/{accession}_2.fastq",
    wrapper:
        "0.72.0/bio/sra-tools/fasterq-dump" 

rule get_small_demo_dataset: 
    input: 
        HTTP.remote("http://genomedata.org/rnaseq-tutorial/HBR_UHR_ERCC_ds_5pc.tar")
    params:
        download_dir="../data/raw/reads"
    output:
        "../data/raw/reads/demo_R1.fq",
        "../data/raw/reads/demo_R2.fq"
    shell:
        "tar -xvf {input} -C {params.download_dir} "
        " && zcat {params.download_dir}/UHR_Rep1_ERCC-Mix1_Build37-ErccTranscripts-chr22.read1.fastq.gz > {output[0]} "
        " && zcat {params.download_dir}/UHR_Rep1_ERCC-Mix1_Build37-ErccTranscripts-chr22.read2.fastq.gz > {output[1]} "
