from pathlib import Path
from tempfile import TemporaryDirectory
from snakemake.shell import shell

with TemporaryDirectory(dir=Path.cwd()) as tmpdir:
    shell("spades.py --rna "
          " --pe1-1 {snakemake.input[0]} "
          " --pe1-2 {snakemake.input[1]} "
          " -o {tmpdir}/ ")
    (Path(tmpdir)/"soft_filtered_transcripts.fasta").rename(snakemake.output)
