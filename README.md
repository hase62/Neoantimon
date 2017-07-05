 #Updated on 11, June. 2017. 
==============================
##1. Preparation
------------------------------
**Set netMHCpan:**
At first, download netMHCpan3.0 from http://www.cbs.dtu.dk/cgi-bin/nph-sw_request?netMHCpan. 
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

**Set netMHCIIpan:**
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
tar -xvf data.tar.gz
```

**Install samtools:**
```
wget http://sourceforge.net/projects/samtools/files/samtools/1.3/samtools-1.3.tar.bz2
tar jxf samtools-1.3.tar.bz2
cd samtools-1.3
./configure
make
make install
cd ..
```

**Install bcftools:**
```
wget https://github.com/samtools/bcftools/releases/download/1.3.1/bcftools-1.3.1.tar.bz2
tar jxf bcftools-1.3.1.tar.bz2
cd bcftools-1.3.1
make
sudo make install
cd ..
```

**Download SampleFiles and CCFP.jar:**
```
wget https://github.com/hase62/Neoantimon/raw/master/lib/ccfp.jar
wget https://github.com/hase62/Neoantimon/raw/master/lib/data.txt.zip
unzip data.txt.zip
```

**Choose Either One of GRCh38, hg38, GRCh37 or hg19...for your sequencing data**
**Download refMrna Files(GRCh38/hg38):**
```
#refMrna Files
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/refMrna.fa.gz
gunzip refMrna.fa.gz
grep ">" refMrna.fa | sed -e "s/>//g" | cut -d' ' -f1 > refMrna.merge.cut1.fa
grep ">" refMrna.fa | sed -e "s/>//g" | cut -d' ' -f2 > refMrna.merge.cut2.fa
grep -v ">" refMrna.fa | sed -e "s/>//g" | cut -d' ' -f2 > refMrna.merge.cut3.fa
paste refMrna.merge.cut1.fa refMrna.merge.cut2.fa refMrna.merge.cut3.fa > refMrna.merge.fa
```

**Download refMrna Files(GRCh37/hg19):**
```
#refMrna Files
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/refMrna.fa.gz
gunzip refMrna.fa.gz
grep ">" refMrna.fa | sed -e "s/>//g" | cut -d' ' -f1 > refMrna.merge.cut1.fa
grep ">" refMrna.fa | sed -e "s/>//g" | cut -d' ' -f2 > refMrna.merge.cut2.fa
grep -v ">" refMrna.fa | sed -e "s/>//g" | cut -d' ' -f2 > refMrna.merge.cut3.fa
paste refMrna.merge.cut1.fa refMrna.merge.cut2.fa refMrna.merge.cut3.fa > refMrna.merge.fa
```

**Download refFlat Files(GRCh38/hg38)**
```
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/refFlat.txt.gz
gunzip refFlat.txt.gz
cut -f2 refFlat.txt > refFlat.cut.txt
```

**Download refFlat Files(GRCh37/hg19)**
```
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/refFlat.txt.gz
gunzip refFlat.txt.gz
cut -f2 refFlat.txt > refFlat.cut.txt
```

**Download human refSeq (GRCh38):**
```
wget ftp://ftp.ensembl.org/pub/release-87/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.toplevel.fa.gz
mv Homo_sapiens.GRCh38.dna.toplevel.fa.gz GRCh38.fa.gz
gunzip GRCh38.fa.gz
samtools faidx GRCh38.fa
```

**Download human refSeq (hg38):**
```
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
gunzip hg38.fa.gz
samtools faidx hg38.fa
```

**Download human refSeq (GRCh37):**
```
wget ftp://ftp.ensembl.org/pub/release-75//fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.75.dna.toplevel.fa.gz
mv Homo_sapiens.GRCh37.75.dna.toplevel.fa.gz GRCh37.fa.gz
gunzip GRCh37.fa.gz
samtools faidx GRCh37.fa
```

##2. Use on R
------------------------------
```
install.packages("devtools")
library(devtools)
install_github('hase62/Neoantimon')
library(Neoantimon)
```

##3. Data Format
------------------------------
```r
library(Neoantimon)
data("hla_table")
print(hla_table, row.names = FALSE)
```

```
##     Name      A1      A2      B1      B2      C1      C2
##   sample A*02:01 A*32:01 B*15:17 B*51:01 C*07:01 C*15:02
##  sample2 A*02:01 A*32:01 B*15:17 B*51:01 C*07:01 C*15:02
```


```r
data("hla_table2")
print(hla_table2, row.names = FALSE)
```

```
##     Name      DPA11      DPA12      DPB11      DPB12      DQA11      DQA12
##   sample DPA1*01:03 DPA1*02:01 DPB1*02:01 DPB1*09:01 DQA1*01:02 DQA1*05:05
##  sample2 DPA1*01:03 DPA1*02:01 DPB1*02:01 DPB1*09:01 DQA1*01:02 DQA1*05:05
##       DQB11      DQB12      DRB11      DRB12
##  DQB1*03:01 DQB1*06:04 DRB1*11:04 DRB1*13:02
##  DQB1*03:01 DQB1*06:04 DRB1*11:04 DRB1*13:02
```


```r
data("sample")
print(sample, row.names = FALSE)
```

```
##  Chr     Start       End Ref Alt Func.refGene Gene.refGene
##    1  47399872  47399872   A   G       exonic      CYP4A11
##    1 116941338 116941338   T   C       exonic       ATP1A1
##    4  24556416  24556416   T   C       exonic        DHX15
##    4  70156404  70156404   -   T       exonic      UGT2B28
##    6  75899298  75899298   T   -       exonic      COL12A1
##    9  89561162  89561162   C   T       exonic         GAS1
##   12  15132141  15132141   G   T       exonic        PDE6H
##   12  20876048  20876048   -   G       exonic      SLCO1C1
##  GeneDetail.refGene ExonicFunc.refGene
##       nonsynonymous                SNV
##          synonymous                SNV
##       nonsynonymous                SNV
##          frameshift          insertion
##          frameshift           deletion
##       nonsynonymous                SNV
##       nonsynonymous                SNV
##          frameshift          insertion
##                                                                                                                                                                                           AAChange.refGene
##                                                                                                                                                                   CYP4A11:NM_000778:exon8:c.T1064C:p.L355P
##                                                                           ATP1A1:NM_000701:exon16:c.T2220C:p.D740D,ATP1A1:NM_001160233:exon16:c.T2220C:p.D740D,ATP1A1:NM_001160234:exon16:c.T2127C:p.D709D
##                                                                                                                                                                     DHX15:NM_001358:exon5:c.A1012G:p.T338A
##                                                                                                                                                                UGT2B28:NM_053039:exon5:c.1186dupT:p.L395fs
##                                                                                                                                                                 COL12A1:NM_004370:exon6:c.628delA:p.I210fs
##                                                                                                                                                                       GAS1:NM_002048:exon1:c.G533A:p.R178H
##                                                                                                                                                                       PDE6H:NM_006205:exon3:c.G163T:p.G55W
##  SLCO1C1:NM_001145944:exon7:c.692_693insG:p.L231fs,SLCO1C1:NM_001145945:exon9:c.899_900insG:p.L300fs,SLCO1C1:NM_017435:exon9:c.1046_1047insG:p.L349fs,SLCO1C1:NM_001145946:exon10:c.1046_1047insG:p.L349fs
##  cytoBand depth_tumor variantNum_tumor depth_normal variantNum_normal
##      1p33          64               28           49                 0
##    1p13.1         100               39          111                 0
##    4p15.2         143               47          151                 0
##    4q13.2          43               15           41                 0
##      6q13         122               38           73                 0
##   9q21.33          20                5           26                 0
##   12p12.3          81               10           47                 0
##   12p12.2          97               11           57                 0
##  bases_tumor bases_normal A_C_G_T_tumor A_C_G_T_normal misRate_tumor
##   46,20,18,8    28,0,21,0     36,0,28,0       49,0,0,0         0.438
##  74,29,26,10    82,0,29,0     0,39,0,61      0,0,0,111         0.390
##  112,39,31,8   129,0,22,0     0,47,0,96      0,0,0,151         0.329
##   28,11,15,4    29,0,12,0           ---            ---         0.349
##   98,31,24,7    63,0,10,0           ---            ---         0.311
##     6,2,14,3    10,0,16,0      0,15,0,5       0,26,0,0         0.250
##    63,8,18,2    35,0,12,0     0,0,71,10       0,0,47,0         0.123
##   82,10,15,1     50,0,7,0           ---            ---         0.113
##  strandRatio_tumor misRate_normal strandRatio_normal P.value.fisher.
##              0.714              0                ---           8.321
##              0.744              0                ---          14.755
##              0.830              0                ---          16.734
##              0.733              0                ---           4.838
##              0.816              0                ---           8.696
##              0.400              0                ---           1.947
##              0.800              0                ---           1.877
##              0.909              0                ---           2.139
##  readPairNum_tumor variantPairNum_tumor otherPairNum_tumor
##                 42                   29                  0
##                 61                   39                  1
##                 98                   50                  0
##                 27                   22                  1
##                 72                   32                  4
##                 15                    5                  0
##                 71                   10                  0
##                 76                   11                  0
##  readPairNum_normal variantPairNum_normal otherPairNum_normal
##                  56                     0                   0
##                 111                     0                   0
##                 152                     0                   1
##                  41                     0                   2
##                  66                     0                   2
##                  27                     0                   0
##                  47                     0                   0
##                  51                     0                   2
##  P.value.fisher_realignment. indel_mismatch_rate indel_mismatch_rate.1
##                        8.617                   0                     0
##                       14.755                   0                     0
##                       17.543                   0                     0
##                        6.926                   0                     0
##                        7.656                   0                     0
##                        1.995                   0                     0
##                        1.877                   0                     0
##                        2.152                   0                     0
##  bp_mismatch_count distance_from_breakpoint simple_repeat_pos
##                  0                        0               ---
##                  0                        0               ---
##                  0                        0               ---
##                  0                        0               ---
##                  0                        0               ---
##                  0                        0               ---
##                  0                        0               ---
##                  1                       16               ---
##  simple_repeat_seq P.value.EBCall.
##                ---          60.000
##                ---          60.000
##                ---          60.000
##                ---          10.783
##                ---          13.076
##                ---           5.208
##                ---           7.481
##                ---           6.179
```


```r
data("RNAseq")
print(RNAseq, row.names = FALSE)
```

```
##  gene_short_name                               locus expression
##              7SK    HSCHR6_MHC_MCF:30910595-30910898 0.00000000
##              7SK    HSCHR6_MHC_QBL:30821624-30821927 0.00000000
##              7SK                 X:12632121-12632316 0.00000000
##             A1BG                19:58856543-58864865 0.19541613
##         A1BG-AS1                19:58859116-58866549 5.30484229
##             A1CF                10:52559168-52645435 4.78487242
##              A2M                  12:9220259-9268825 0.75619084
##          A2M-AS1                  12:9217772-9220651 4.15985770
##            A2ML1                  12:8975067-9039597 7.50201326
##        A2ML1-AS1                  12:8928814-8983543 0.00000000
##        A2ML1-AS2                  12:8972411-8973309 7.02685468
##            A2MP1                  12:9381128-9428413 4.36764217
##          A3GALT2                 1:33772366-33786699 0.00000000
##           A4GALT                22:43088126-43117304 2.27225123
##            A4GNT               3:137842559-137851229 0.00000000
##             AAAS                12:53701239-53718648 1.45197987
##             AACS              12:125549924-125627873 7.34085879
##           AACSP1               5:178191861-178245436 4.07342571
##            AADAC               3:151531824-151546276 0.61589770
##          AADACL2               3:151451703-151479127 0.01359431
##          AADACL2 HSCHR3_1_CTG2_1:151462241-151489665 0.00000000
```


```r
data("CopyNum")
print(CopyNum, row.names = FALSE)
```

```
##  Chromosome Position      Log.R segmented.LogR BAF segmented.BAF
##           1   564621  0.6071447    -0.09862298   1            NA
##           1   799463  0.1519967    -0.09862298   1            NA
##           1  1017216  0.8146911    -0.09862298   0            NA
##           1  1158277 -1.9594627    -0.09862298   0            NA
##           1  1242215  0.2927962    -0.09862298   0            NA
##           1  1462766 -0.2234090    -0.09862298   1            NA
##  Copy.number Minor.allele Raw.copy.number
##            2            1       4.3540752
##            2            1       2.1467339
##            3            1       5.8658499
##            1            0      -0.5035897
##            2            1       2.6999875
##            2            1       1.0726687
```

##4. Sample Codes
------------------------------
##
#./lib/
ccfp.jar
GRCh37.fa
GRCh37.fa.fai
netMHCIIpan-3.1
netMHCpan-3.0
refFlat.cut.txt
refFlat.txt
refMrna.fa
refMrna.merge.cut1.fa
refMrna.merge.cut2.fa
refMrna.merge.cut3.fa
refMrna.merge.fa
#
#./data/
CopyNum.txt
hla_table.txt
hla_table2.txt
RNAseq.txt
sample.txt
##

Calculate A List of Neoantigens on SNVs for HLA Class I. 
```
MainSNVClass1(hmdir = getwd(),
              input_file = "data/sample.txt",
              job_ID = "NO_JOB_ID",
              Chr_Column = 1,
              Mutation_Start_Column = 2,
              Mutation_End_Column = 3,
              Mutation_Ref_Column = 4,
              Mutation_Alt_Column = 5,
              NM_ID_Column = 10,
              Depth_Tumor_Column = 12,
              Depth_Normal_Column = 14,
              file_name_in_HLA_table = "sample",
              HLA_file = "data/hla_table.txt",
              RNAseq_file = "data/RNAseq.txt",
              CNV="data/Copy.txt",
              Purity = 0.8,
              ccfp_dir = "lib/ccfp.jar",
              netMHCpan_dir = "lib/netMHCpan-3.0/netMHCpan",
              refDNA = "lib/GRCh37.fa")
```

Marge Results for HLA Class I. 
```
MainMergeClass1(hmdir = getwd(),
                input_dir = "data",
                input_file_prefix = "sample",
                Tumor_RNA_BASED_ON_DNA = TRUE)
```

Calculate A List of Neoantigens for HLA Class I. 
```
MainSNVClass2(hmdir = getwd(),
              input_file = "data/sample.txt",
              job_ID = "NO_JOB_ID",
              Chr_Column = 1,
              Mutation_Start_Column = 2,
              Mutation_End_Column = 3,
              Mutation_Ref_Column = 4,
              Mutation_Alt_Column = 5,
              NM_ID_Column = 10,
              Depth_Tumor_Column = 12,
              Depth_Normal_Column = 14,
              file_name_in_HLA_table = "sample",
              HLA_file = "data/hla_table2.txt",
              RNAseq_file = "data/RNAseq.txt",
              CNV="data/Copy.txt",
              Purity = 0.8,
              ccfp_dir = "lib/ccfp.jar",
              netMHCpan_dir = "lib/netMHCIIpan-3.1/netMHCIIpan",
              refDNA = "lib/GRCh37.fa")
```

Marge Results for HLA Class II. 
```
MainMergeClass2(hmdir = getwd(),
                input_dir = "data",
                input_file_prefix = "sample",
                Tumor_RNA_BASED_ON_DNA = TRUE)
```

Calculate A List of Neoantigens on Indels for HLA Class I. 
```
MainINDELClass1(hmdir = getwd(),
                input_file = "data/sample.txt",
                job_ID = "NO_JOB_ID",
                Chr_Column = 1,
                Mutation_Start_Column = 2,
                Mutation_End_Column = 3,
                Mutation_Ref_Column = 4,
                Mutation_Alt_Column = 5,
                NM_ID_Column = 10,
                Depth_Tumor_Column = 12,
                Depth_Normal_Column = 14,
                file_name_in_HLA_table = "sample",
                HLA_file = "data/hla_table.txt",
                RNAseq_file = "data/RNAseq.txt",
                CNV="data/Copy.txt",
                Purity = 0.8,
                ccfp_dir = "lib/ccfp.jar",
                netMHCpan_dir = "lib/netMHCpan-3.0/netMHCpan",
                refDNA = "lib/GRCh37.fa")
```

Calculate A List of Neoantigens on Indels for HLA Class II. 
```
MainINDELClass2(hmdir = getwd(),
                input_file = "data/sample.txt",
                job_ID = "NO_JOB_ID",
                Chr_Column = 1,
                Mutation_Start_Column = 2,
                Mutation_End_Column = 3,
                Mutation_Ref_Column = 4,
                Mutation_Alt_Column = 5,
                NM_ID_Column = 10,
                Depth_Tumor_Column = 12,
                Depth_Normal_Column = 14,
                file_name_in_HLA_table = "sample",
                HLA_file = "data/hla_table2.txt",
                RNAseq_file = "data/RNAseq.txt",
                CNV="data/Copy.txt",
                Purity = 0.8,
                ccfp_dir = "lib/ccfp.jar",
                netMHCpan_dir = "lib/netMHCIIpan-3.1/netMHCIIpan",
                refDNA = "lib/GRCh37.fa")
```




