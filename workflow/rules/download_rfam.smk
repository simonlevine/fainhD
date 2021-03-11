from snakemake.remote.FTP import RemoteProvider as FtpRemoteProvider
FTP = FtpRemoteProvider()

rule download_rfam_db:
    input:
        FTP.remote("ftp://ftp.ebi.ac.uk/pub/databases/Rfam/14.3/Rfam.cm.gz")
    output:
        "data/raw/rfam/Rfam.cm"
    conda:
        "envs/infernal.yaml"
    shell:
        "zcat {input} > {output} && cmpress {output}"