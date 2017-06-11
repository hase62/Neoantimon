#!/bin/bash
#$ -S /bin/bash
#$ -cwd

#Make refMrna Files
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/refMrna.fa.gz
gunzip refMrna.fa.gz
grep ">" refMrna.fa | sed -e "s/>//g" | cut -d' ' -f1 > refMrna.merge.cut1.fa
grep ">" refMrna.fa | sed -e "s/>//g" | cut -d' ' -f2 > refMrna.merge.cut2.fa
grep -v ">" refMrna.fa | sed -e "s/>//g" | cut -d' ' -f2 > refMrna.merge.cut3.fa
paste refMrna.merge.cut1.fa refMrna.merge.cut2.fa refMrna.merge.cut3.fa > refMrna.merge.fa

#Make refFlat Files
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/refFlat.txt.gz
gunzip refFlat.txt.gz
cut -f2 refFlat.txt > refFlat.cut.txt
