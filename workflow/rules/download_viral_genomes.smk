from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider
HTTP = HTTPRemoteProvider()

rule download_virus_genomes:
    input:
        HTTP.remote("https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.1.1.genomic.fna.gz"),
    output:
        "data/raw/viral_genomes.fasta"
    shell:
        "gzip -d -c {input} >| {output}"
