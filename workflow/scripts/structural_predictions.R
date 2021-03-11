# if (!requireNamespace("BiocManager", quietly = TRUE))
    # install.packages("BiocManager", repos='http://cran.us.r-project.org')
# if (!requireNamespace("tidyverse", quietly = TRUE))
	# install.packages("tidyverse", repos='http://cran.us.r-project.org')
# The following initializes usage of Bioc devel
# BiocManager::install(version='devel')
# BiocManager::install("rfaRm")
# BiocManager::install("Biostrings")

library(tidyverse)
library(rfaRm)
library(Biostrings)

# init our data dumps
RNAfamilies <- tibble(originalSequence = "", rfamAccession= "",rfamID="", eValue="", alignmentMatch="", alignmentSecondaryStructure="", rfamFamilySummary="")

ORFS <- vector()
convertedORFS <- vector()


# parse the fasta and peel off each ORF, store as list of ORF sequences
# data/processed/{sample}_orf.fasta
ORFS <- readDNAStringSet(snakemake@input[[1]])

for (ORF in length(ORFs)) {
	# convert ORF DNA to mRNA
	append(convertedORFS, ReverseComplement(ORFs[ORF]))
	}

# search each mRNA sequence for RNA structure
for (sequence in length(convertedORFS)){
	SearchHits <- rfamSequenceSearch(convertedORFS[sequence])
	for (hit in length(SearchHits)){
		RNAfamilies <-add_row(
		convertedORFS[sequence],
		SearchHits[[hit]]$rfamAccession,
		SearchHits[[hit]]$rfamID,
		SearchHits[[hit]]$eValue,
		SearchHits[[hit]]$alignmentMatch,
		SearchHits[[hit]]$alignmentSecondaryStructure,
		rfamFamilySummary(SearchHits[[hit]]$rfamAccession)
	) }}

# output RNAfamilies as tsv
write_csv(RNAfamilies, snakemake@output[[1]]) 




