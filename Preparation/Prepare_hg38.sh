#!/bin/bash
#$ -S /bin/bash
#$ -cwd

#Download Ref Data
mkdir lib
cd lib
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz

#Unpack
gunzip hg38.fa.gz

#Make Index
samtools faidx hg38.fa
