from pathlib import Path
from tempfile import TemporaryDirectory
from snakemake.shell import shell

with TemporaryDirectory(dir=Path.cwd()) as tmpdir:
    shell("cmscan --rfam --cut_ga --nohmmonly"
          "--tblout {output[0]} " 
          "{input.rfam_database} "
          "{input.fasta} "
          "> {output[1]}")
    (Path(tmpdir)/"soft_filtered_transcripts.fasta") \
        .rename(snakemake.output[0])



#cmscan --rfam --cut_ga --nohmmonly --tblout mrum-genome.tblout Rfam.cm contigs.fasta > myfavscan.cmscan             
 