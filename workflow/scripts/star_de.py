from pathlib import Path
from tempfile import TemporaryDirectory
from snakemake.shell import shell

with TemporaryDirectory() as tmpdir:
    shell("STAR --runThreadN {snakemake.threads} "
          "--genomeDir {snakemake.input.reference_genome_dir} "
          "--readFilesIn {snakemake.input.fq1} {snakemake.input.fq2} "
          "--outFileNamePrefix {tmpdir} "
          "--outReadsUnmapped Fastx ")
    for outpath, fname in [(snakemake.output[0], "Unmapped.out.mate1"),
                           (snakemake.output[1], "Unmapped.out.mate2")]:
        (Path(tmpdir)/fname).rename(outpath)
