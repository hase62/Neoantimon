#!/bin/bash
#$ -S /bin/bash
#$ -cwd

#Download Ref Data
mkdir lib
cd lib
wget ftp://ftp.ensembl.org/pub/release-87/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.toplevel.fa.gz
mv Homo_sapiens.GRCh38.dna.toplevel.fa.gz GRCh38.fa.gz

#Unpack
gunzip GRCh38.fa.gz

#Make Index
samtools faidx GRCh38.fa.gz
