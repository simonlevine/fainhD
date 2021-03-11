import pandas as pd
from Bio import SeqIO

viral_contigs = pd.read_csv(snakemake.input[0], header=None, delimiter="\t") \
    .iloc[:,0].unique()
viral_records = [
    contig for contig in SeqIO.parse(snakemake.input[1], "fasta")
    if contig.id in viral_contigs]
with open(snakemake.output[0], "wt") as fout: 
    SeqIO.write(viral_records, fout, "fasta")