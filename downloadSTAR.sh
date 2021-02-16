#!/bin/sh

url='http://labshare.cshl.edu/shares/gingeraslab/www-data/dobin/STAR/STARgenomes/Human/GRCh38_Ensembl99_sparseD3_sjdbOverhang99/'

for f in {'chrLength.txt' 'chrName.txt' 'chrStart.txt' 'Genome' 'genomeParameters.txt' 'SA' 'SAindex'}; do
echo $f
wget 'http://labshare.cshl.edu/shares/gingeraslab/www-data/dobin/STAR/STARgenomes/Human/GRCh38_Ensembl99_sparseD3_sjdbOverhang99/'$f -O 'input/reference/'$f
done