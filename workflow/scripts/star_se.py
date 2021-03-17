from pathlib import Path
from tempfile import TemporaryDirectory
from snakemake.shell import shell

read_proc = "zcat" if snakemake.input.fq.endswith(".gz") else ""

with TemporaryDirectory(dir=Path.cwd()) as tmpdir:
    shell("STAR --runThreadN {snakemake.threads} "
          "--genomeDir {snakemake.input.reference_genome_dir} "
          "--readFilesIn {snakemake.input.fq} "
          "{read_proc} "
          "--outFileNamePrefix {tmpdir}/ "
          "--outReadsUnmapped Fastx ")
    (Path(tmpdir)/"Unmapped.out.mate1").rename(snakemake.output[0])
