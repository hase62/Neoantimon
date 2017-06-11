#!/bin/bash
#$ -S /bin/bash
#$ -cwd

#Install bcftools
wget https://github.com/samtools/bcftools/releases/download/1.3.1/bcftools-1.3.1.tar.bz2
tar jxf bcftools-1.3.1.tar.bz2
cd bcftools-1.3.1
make
sudo make install
cd ..
