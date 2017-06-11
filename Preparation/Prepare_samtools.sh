#!/bin/bash
#$ -S /bin/bash
#$ -cwd

#Install samtools
wget http://sourceforge.net/projects/samtools/files/samtools/1.3/samtools-1.3.tar.bz2
tar jxf samtools-1.3.tar.bz2
cd samtools-1.3
./configure
make
make install
cd ..