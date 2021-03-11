import re
from functools import reduce
import pandas as pd
from Bio import SeqIO


def main():
    dfs = [ingest_structural_predictions(snakemake.input[0]),
           ingest_blast_results(snakemake.input[1]),
           ingest_contigs(snakemake.input[2]),
           ingest_orfs(snakemake.input[3])]
    df = reduce(lambda a, b: pd.merge(a, b, on="query_name", how="outer"), dfs)
    df.to_json(snakemake.output[0], lines=True, orient="records")


def ingest_structural_predictions(fp) -> pd.DataFrame:
    records = []
    with open(fp) as f:
        for i, line in enumerate(f):
            if not line.startswith("#"):
                line = [c.replace("\n", "") for c in line.split(" ") if c]
                target_name, accession, query_name, accession2, mdl, mdl_from, \
                mdl_to, seq_from, sqe_to, strand, trunc, pass_, gc, bias, score, \
                evalue, inc, *desc = line
                records.append([query_name, target_name, evalue, " ".join(desc)])
    df = pd.DataFrame(records)
    df.columns = ["query_name", "target_name", "rfam_e_value", "structural_description"]
    return df


def ingest_blast_results(fp) -> pd.DataFrame:
    df_blast = pd.read_csv(fp) \
        .rename(columns={"qseqid": "query_name", "evalue": "blast_e_value"})
    df_blast_linewise = df_blast.groupby("query_name").apply(
        lambda df: [
            {"sseqid": sseqid, "evalue": evalue}
            for (sseqid, evalue) in zip(df.sseqid, df.blast_e_value)
        ]
    ).to_frame().rename(columns={0: "blast_results"}).reset_index()
    return df_blast_linewise


def ingest_contigs(fp) -> pd.DataFrame:
    records = []
    for record in SeqIO.parse(fp, "fasta"):
        records.append({"query_name": record.id, "contig_seq": str(record.seq)})
    contig_df = pd.DataFrame(records)
    return contig_df


def ingest_orfs(fp) -> pd.DataFrame:
    records = []
    for record in SeqIO.parse(fp, "fasta"):
        query_name, = re.match(r"^(.*)_ORF\.\d", record.id).groups()
        records.append({"query_name": query_name, "orfs": str(record.seq)})
    orf_df = pd.DataFrame(records).groupby("query_name").orfs.apply(list).to_frame()
    return orf_df


if __name__ == "__main__":
    main()