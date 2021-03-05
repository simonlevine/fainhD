from snakemake.remote.HTTP import RemoteProvider
HTTP = RemoteProvider()

rule download_from_SRA:
    output:
        "data/raw/reads/{accession}_1.fastq",
        "data/raw/reads/{accession}_2.fastq",
    wrapper:
        "0.72.0/bio/sra-tools/fasterq-dump" 