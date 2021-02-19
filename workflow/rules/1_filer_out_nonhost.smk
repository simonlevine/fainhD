from snakemake.remote.HTTP import RemoteProvider

HTTP = RemoteProvider()
reference_genome_url_prefix = "http://labshare.cshl.edu/shares/gingeraslab/www-data/dobin/STAR/STARgenomes/Human/GRCh38_Ensembl99_sparseD3_sjdbOverhang99"

rule download_genome:
    input:
        [HTTP.remote(f"{reference_genome_url_prefix}/{f}", keep_local=True)
         for f in ['chrLength.txt', 'chrName.txt', 'chrStart.txt', 'Genome',
                   'genomeParameters.txt', 'SA', 'SAindex', "sjdbInfo.txt"]]
    output:
        directory("../data/raw/reference_genome")
    run:
        shell("mkdir -p {output}")
        for f in input:
            shell("mv {f} {output}")

rule star_double_ended:
    input:
        fq1 = "../data/raw/reads/{sample}_R1.fq",
        fq2 = "../data/raw/reads/{sample}_R2.fq",
        reference_genome_dir=rules.download_genome.output[0]
    output:
        "../data/interim/{sample}_aligned_to_human_genome.sam",
    conda:
        "../envs/star.yaml"
    threads: 8
    shell:
        "mkdir -p ../data/interim/{wildcards.sample}_alignment_workingdir/ "
        "&& STAR"
        " --runThreadN {threads}"
        " --genomeDir {input.reference_genome_dir}"
        " --readFilesIn {input.fq1} {input.fq2}"
        " --outFileNamePrefix ../data/interim/{wildcards.sample}_alignment_workingdir/ "
        " --outSAMunmapped Within"
        "&& mv ../data/interim/{wildcards.sample}_alignment_workingdir/Aligned.out.sam {output} "
        "&& rm -rf ../data/interim/{wildcards.sample}_alignment_workingdir/"

rule extract_unmapped_reads:
   input:
        "../data/interim/{sample}_aligned_to_human_genome.sam"
   output:
        "../data/interim/{sample}_nonhost.sam"
   params:
       # '4' is the flag for unmapped reads
       # see: http://www.htslib.org/doc/samtools.html
       "-f 4"
   wrapper:
       "0.72.0/bio/samtools/view" 

rule convert_sam_to_fastq:
    input:
        "../data/interim/{sample}_nonhost.sam"
    output:
        "../data/interim/{sample}_nonhost.fq"
    run:
        with open(input[0]) as fin, open(output[0], "w") as fout:
            records = fin.readlines()
            for i, record in enumerate(records):
                if i % 10000:
                    print(f"converting sam to fastq: {i/len(records):.2%}...", flush=True)
                qname, flag, rame, pos, mapq, cigar, rnext, pnetx, tlen, seq, qual, *_ = record.split("\t")
                fout.write(f"@{qname}\n{seq}\n+\n{qual}\n")
