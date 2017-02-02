#!/bin/bash
#$ -S /bin/bash
#$ -cwd

#
mkdir ./../lib_int
cd ./../lib_int
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/refFlat.txt.gz
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/refMrna.fa.gz
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
wget ftp://ftp.ensembl.org/pub/release-87/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.toplevel.fa.gz
mv Homo_sapiens.GRCh38.dna.toplevel.fa GRCh38.fa

#
gunzip refFlat.txt.gz
gunzip refMrna.fa.gz
gunzip hg38.fa.gz
gunzip Homo_sapiens.GRCh38.dna.toplevel.fa.gz

#
samtools faidx hg38.fa
samtools Homo_sapiens.GRCh38.dna.toplevel.fa

#
if [ -e netMHCpan-3.0a.Linux.tar ]; then
 tar zxvf netMHCIIpan-3.1a.Linux.tar
else
 mv ./../main/netMHCpan-3.0a.Linux.tar ./
 tar zxvf netMHCpan-3.0a.Linux.tar 
fi

#
mkdir netMHCpan-3.0/tmp
sed -i -e "s/\/usr\/cbs\/packages\/netMHCpan\/3.0\/netMHCpan-3.0/${0}/g" netMHCpan-3.0/netMHCpan
sed -i -e "s/#setenv/setenv/g" netMHCpan-3.0/netMHCpan

#
if [ -e netMHCIIpan-3.1a.Linux.tar ]; then
 tar zxvf netMHCIIpan-3.1a.Linux.tar
else
 mv ./../main/netMHCIIpan-3.1a.Linux.tar ./
 tar zxvf netMHCIIpan-3.1a.Linux.tar 
fi

#
mkdir netMHCIIpan-3.1/tmp
sed -i -e "s/\/usr\/cbs\/packages\/netMHCIIpan\/3.1\/netMHCIIpan-3.1/${0}/g" netMHCIIpan-3.1/netMHCIIpan
sed -e "15i setenv  TMPDIR  \$\{NMHOME\}/tmp" netMHCIIpan-3.1/netMHCIIpan

#
grep ">" refMrna.fa | sed -e "s/>//g" | cut -d' ' -f1 > refMrna.merge.cut1.fa_
grep ">" refMrna.fa | sed -e "s/>//g" | cut -d' ' -f2 > refMrna.merge.cut2.fa_
grep -v ">" refMrna.fa | sed -e "s/>//g" | cut -d' ' -f2 > refMrna.merge.cut3.fa_
paste refMrna.merge.cut1.fa refMrna.merge.cut2.fa refMrna.merge.cut3.fa > refMrna.merge.fa

#
cd ./../main