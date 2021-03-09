from pathlib import Path
from tempfile import TemporaryDirectory
from snakemake.shell import shell

read_proc = ""  # plain text fastqs do not need a process to interface with STAR
compressed_fastq_as_input = snakemake.input.fq1.endswith(".gz")
if compressed_fastq_as_input:
      assert snakemake.input.fq2.endswith(".gz"), \
            "paired read files may be compressed or uncompressed, but must match!"
      read_proc = "--readFilesIn zcat"
else:
      assert not snakemake.input.fq2.endswith(".gz"), 
            "paired read files may be compressed or uncompressed, but must match!"

with TemporaryDirectory(dir=Path.cwd()) as tmpdir:
    shell("STAR --runThreadN {snakemake.threads} "
          "--genomeDir {snakemake.input.reference_genome_dir} "
          "--readFilesIn {snakemake.input.fq1} {snakemake.input.fq2} "
          "{read_proc} "
          "--outFileNamePrefix {tmpdir}/ "
          "--outReadsUnmapped Fastx ")
    for outpath, fname in [(snakemake.output[0], "Unmapped.out.mate1"),
                           (snakemake.output[1], "Unmapped.out.mate2")]:
        (Path(tmpdir)/fname).rename(outpath)
