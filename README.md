#Updated on 22th, Feb. 2017. 
==============================
##Preparation
------------------------------
**Set netMHCpan**
Download netMHCpan3.0 from http://www.cbs.dtu.dk/cgi-bin/nph-sw_request?netMHCpan. 
```
#Set netMHCpan3.0 bash

#!/bin/bash
#$ -S /bin/bash
#$ -cwd

if [ -e netMHCpan-3.0a.Linux.tar.gz ]; then
 tar zxvf netMHCpan-3.0a.Linux.tar.gz
fi
if [ -e netMHCpan-3.0a.Darwin.tar.gz ]; then
 tar zxvf netMHCpan-3.0a.Darwin.tar.gz
fi
cd netMHCpan-3.0
mkdir tmp
cdir=`pwd`
sed -i -e "s:/usr/cbs/packages/netMHCpan/3.0/netMHCpan-3.0:${cdir}:g" netMHCpan
sed -i -e "s:#setenv:setenv:g" netMHCpan
wget http://www.cbs.dtu.dk/services/NetMHCpan-3.0/data.tar.gz
tar -xvf data.tar.gz
```

**Set netMHCIIpan**
Download netMHCIIpan 3.1 from http://www.cbs.dtu.dk/cgi-bin/nph-sw_request?netMHCIIpan. 
```
#Set netMHCpan3.0 bash

#!/bin/bash
#$ -S /bin/bash
#$ -cwd

if [ -e netMHCIIpan-3.1a.Linux.tar.gz ]; then
 tar zxvf netMHCIIpan-3.1a.Linux.tar.gz
elif [ -e netMHCIIpan-3.1a.Darwin.tar.gz ]
 tar zxvf netMHCIIpan-3.1a.Darwin.tar.gz
fi
cd netMHCIIpan-3.1
mkdir tmp
cdir=`pwd`
sed -i -e "s:/usr/cbs/packages/netMHCIIpan/3.1/netMHCIIpan-3.1:${cdir}:g" netMHCIIpan
sed -i -e "15s:^:setenv  TMPDIR \$\{NMHOME\}/tmp:" netMHCIIpan
wget http://www.cbs.dtu.dk/services/NetMHCIIpan-3.1/data.tar.gz
#gunzip -c data.tar.gz | tar xvf -
tar -xvf data.tar.gz```

**Install samtools**
```
#!/bin/bash
#$ -S /bin/bash
#$ -cwd

wget http://sourceforge.net/projects/samtools/files/samtools/1.3/samtools-1.3.tar.bz2
tar jxf samtools-1.3.tar.bz2
cd samtools-1.3
./configure
make
make install
cd ..
```

**Install bcftools**
```
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
```

**Download refFiles**
```
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
```

**Download human refSeq (hg38)**
```
#!/bin/bash
#$ -S /bin/bash
#$ -cwd

mkdir lib
cd lib
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz

#Unpack
gunzip hg38.fa.gz

#Make Index
samtools faidx hg38.fa
```

**Download human refSeq (GRCh38)**
```
#!/bin/bash
#$ -S /bin/bash
#$ -cwd

mkdir lib
cd lib
wget ftp://ftp.ensembl.org/pub/release-87/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.toplevel.fa.gz
mv Homo_sapiens.GRCh38.dna.toplevel.fa.gz GRCh38.fa.gz
gunzip GRCh38.fa.gz
samtools faidx GRCh38.fa.gz
```

**Download CCFP.jar**
```
#!/bin/bash
#$ -S /bin/bash
#$ -cwd

wget ccfp.jar
```



(2) 
install.packages("devtools")
library(devtools)
install_github('hase62/Neoantimon')
library(Neoantimon)
vignette("SampleCodeForNeoantimon")

(3)See following texts and “SampleRun.R” to learn how to use this library. 
