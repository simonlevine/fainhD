from pathlib import Path
from tempfile import TemporaryDirectory
from snakemake.shell import shell

if snakemake.params["sequencing"] in ["se", "single-ended"]:
    fq1, = snakemake.input
    read_spec = f"--s1 {fq1}"
elif snakemake.params["sequencing"] in ["pe", "pair-ended"]:
    fq1, fq2 = snakemake.input
    read_spec = f"--pe1-1 {fq1} --pe-2 {fq2}"
else:
    raise ValueError("Please specify a read chemistry! ('se' or 'pe')")

with TemporaryDirectory(dir=Path.cwd()) as tmpdir:
    shell("spades.py --rna {read_spec} -o {tmpdir}/ ")
    (Path(tmpdir)/"soft_filtered_transcripts.fasta") \
        .rename(snakemake.output[0])