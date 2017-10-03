## Updated on 3, Oct. 2017. 
## 1. Preparation
**Set netMHCpan:**

1. Download netMHCpan3.0 from http://www.cbs.dtu.dk/cgi-bin/nph-sw_request?netMHCpan. 

2. Make a script (setNetMHCpan.sh) described below or download it from https://github.com/hase62/Neoantimon/raw/master/lib, 
and then run the script as "./setNetMHCpan.sh". 
```
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

1. Download netMHCIIpan 3.1 from http://www.cbs.dtu.dk/cgi-bin/nph-sw_request?netMHCIIpan. 

2. Make a script (setNetMHCIIpan.sh) described below or download it from https://github.com/hase62/Neoantimon/raw/master/lib, 
and then run the script as "./setNetMHCIIpan.sh". 
```
#!/bin/bash
#$ -S /bin/bash
#$ -cwd

if [ -e netMHCIIpan-3.1a.Linux.tar.gz ]; then
 tar zxvf netMHCIIpan-3.1a.Linux.tar.gz
fi
if [ -e netMHCIIpan-3.1a.Darwin.tar.gz ]; then
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

**Install samtools and bcftools:**

Run the following codes or download files from https://github.com/hase62/Neoantimon/raw/master/lib/. 
```
wget http://sourceforge.net/projects/samtools/files/samtools/1.6/samtools-1.6.tar.bz2
tar jxf samtools-1.6.tar.bz2
cd samtools-1.6
./configure
make
make install
cd ..
```

```
wget https://github.com/samtools/bcftools/releases/download/1.6/bcftools-1.6.tar.bz2
tar jxf bcftools-1.6.tar.bz2
cd bcftools-1.6
make
sudo make install
cd ..
```

In addition, if you want to calculate variant allele frequency (VAF), get the old one. 
```
https://sourceforge.net/projects/samtools/files/samtools/0.1.19/samtools-0.1.19.tar.bz2
tar jxf samtools-0.1.19.tar.bz2
cd samtools-0.1.19
make
cd ..
```

**Download SampleFiles and CCFP.jar:**
```
wget https://github.com/hase62/Neoantimon/raw/master/lib/ccfp.jar
wget https://github.com/hase62/Neoantimon/raw/master/lib/data.txt.zip
unzip data.txt.zip
```

**Download refMrna Files (Required, you have to get your corresponding version from GRCh38, hg38, GRCh37 or hg19):**

*Download refMrna Files(GRCh38/hg38)*
```
#refMrna Files
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/refMrna.fa.gz
gunzip refMrna.fa.gz
grep ">" refMrna.fa | sed -e "s/>//g" | cut -d' ' -f1 > refMrna.cut1.fa
grep ">" refMrna.fa | sed -e "s/>//g" | cut -d' ' -f2 > refMrna.cut2.fa
perl -pe 's/[0-9, ,NM_,NR_,\n]//g' refMrna.fa > tmp
perl -pe 's/>/\n/g' tmp > refMrna.cut3.fa
sed -i -e '/^$/d' refMrna.cut3.fa
rm tmp
paste refMrna.cut1.fa refMrna.cut2.fa refMrna.cut3.fa > refMrna.merge.fa
```

*Download refMrna Files(GRCh37/hg19):*
```
#refMrna Files
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/refMrna.fa.gz
gunzip refMrna.fa.gz
grep ">" refMrna.fa | sed -e "s/>//g" | cut -d' ' -f1 > refMrna.cut1.fa
grep ">" refMrna.fa | sed -e "s/>//g" | cut -d' ' -f2 > refMrna.cut2.fa
perl -pe 's/[0-9, ,NM_,NR,\n]//g' refMrna.fa > tmp
perl -pe 's/>/\n/g' tmp > refMrna.cut3.fa
sed -i -e '/^$/d' refMrna.cut3.fa
rm tmp
paste refMrna.cut1.fa refMrna.cut2.fa refMrna.cut3.fa > refMrna.merge.fa
```

**Download refFlat Files (Required, you have to get your corresponding version from GRCh38, hg38, GRCh37 or hg19)**

*Download refFlat Files(GRCh38/hg38)*
```
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/refFlat.txt.gz
gunzip refFlat.txt.gz
cut -f2 refFlat.txt > refFlat.cut.txt
```

*Download refFlat Files(GRCh37/hg19)*
```
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/refFlat.txt.gz
gunzip refFlat.txt.gz
cut -f2 refFlat.txt > refFlat.cut.txt
```

**Download human refSeq (Not required, if you want to calculate variant allele frequency (VAF)):**

*Download human refSeq (GRCh38):*
```
wget ftp://ftp.ensembl.org/pub/release-87/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.toplevel.fa.gz
mv Homo_sapiens.GRCh38.dna.toplevel.fa.gz GRCh38.fa.gz
gunzip GRCh38.fa.gz
samtools faidx GRCh38.fa
```

*Download human refSeq (hg38):*
```
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
gunzip hg38.fa.gz
samtools faidx hg38.fa
```

*Download human refSeq (GRCh37):*
```
wget ftp://ftp.ensembl.org/pub/release-75//fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.75.dna.toplevel.fa.gz
mv Homo_sapiens.GRCh37.75.dna.toplevel.fa.gz GRCh37.fa.gz
gunzip GRCh37.fa.gz
samtools faidx GRCh37.fa
```

## 2. Use on R
```
install.packages("devtools");
library(devtools);
install_github('hase62/Neoantimon');
library(Neoantimon);
```

## 3. Data Format
------------------------------
*A HLA table is required.*
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

*A HLA table is required.*
```r
data("hla_table2")
print(hla_table2, row.names = FALSE)
```

```
##	Name	DPA11	DPA12	DPB11	DPB12	DQA11	DQA12	DQB11	DQB12	DRB11	DRB12
##	sample	DPA1*01:03	DPA1*02:01	DPB1*02:01	DPB1*09:01	DQA1*01:02	DQA1*05:05	DQB1*03:01	DQB1*06:04	DRB1*11:04	DRB1*13:02
##	sample2	DPA1*01:03	DPA1*02:01	DPB1*02:01	DPB1*09:01	DQA1*01:02	DQA1*05:05	DQB1*03:01	DQB1*06:04	DRB1*11:04	DRB1*13:02
```

*An annotated VCF file is required. It must include columns representing "Chromosome Number", "Mutation Start Position", "Mutation End Position", "Mutation Ref", "Mutation Alt", and "NM_ID".*
```r
data("sample")
print(sample, row.names = FALSE)
```

```
##	Chr	Start	End	Ref	Alt	Func.refGene	Gene.refGene	GeneDetail.refGene	ExonicFunc.refGene	AAChange.refGene	cytoBand	depth_tumor	variantNum_tumor	depth_normal	variantNum_normal	bases_tumor	bases_normal	A_C_G_T_tumor	A_C_G_T_normal	misRate_tumor	strandRatio_tumor	misRate_normal	strandRatio_normal	P.value.fisher.	readPairNum_tumor	variantPairNum_tumor	otherPairNum_tumor	readPairNum_normal	variantPairNum_normal	otherPairNum_normal	P.value.fisher_realignment.	indel_mismatch_rate	indel_mismatch_rate.1	bp_mismatch_count	distance_from_breakpoint	simple_repeat_pos	simple_repeat_seq	P.value.EBCall.
##	1	47399872	47399872	A	G	exonic	CYP4A11	nonsynonymous	SNV	CYP4A11:NM_000778:exon8:c.T1064C:p.L355P	1p33	64	28	49	0	46,20,18,8	28,0,21,0	36,0,28,0	49,0,0,0	0.438	0.714	0	---	8.321	42	29	0	56	0	0	8.617	0	0	0	0	---	---	60
##	1	116941338	116941338	T	C	exonic	ATP1A1	synonymous	SNV	ATP1A1:NM_000701:exon16:c.T2220C:p.D740D,ATP1A1:NM_001160233:exon16:c.T2220C:p.D740D,ATP1A1:NM_001160234:exon16:c.T2127C:p.D709D	1p13.1	100	39	111	0	74,29,26,10	82,0,29,0	0,39,0,61	0,0,0,111	0.39	0.744	0	---	14.755	61	39	1	111	0	0	14.755	0	0	0	0	---	---	60
##	4	24556416	24556416	T	C	exonic	DHX15	nonsynonymous	SNV	DHX15:NM_001358:exon5:c.A1012G:p.T338A	4p15.2	143	47	151	0	112,39,31,8	129,0,22,0	0,47,0,96	0,0,0,151	0.329	0.83	0	---	16.734	98	50	0	152	0	1	17.543	0	0	0	0	---	---	60
##	4	70156404	70156404	-	T	exonic	UGT2B28	frameshift	insertion	UGT2B28:NM_053039:exon5:c.1186dupT:p.L395fs	4q13.2	43	15	41	0	28,11,15,4	29,0,12,0	---	---	0.349	0.733	0	---	4.838	27	22	1	41	0	2	6.926	0	0	0	0	---	---	10.783
##	6	75899298	75899298	T	-	exonic	COL12A1	frameshift	deletion	COL12A1:NM_004370:exon6:c.628delA:p.I210fs	6q13	122	38	73	0	98,31,24,7	63,0,10,0	---	---	0.311	0.816	0	---	8.696	72	32	4	66	0	2	7.656	0	0	0	0	---	---	13.076
##	9	89561162	89561162	C	T	exonic	GAS1	nonsynonymous	SNV	GAS1:NM_002048:exon1:c.G533A:p.R178H	9q21.33	20	5	26	0	6,2,14,3	10,0,16,0	0,15,0,5	0,26,0,0	0.25	0.4	0	---	1.947	15	5	0	27	0	0	1.995	0	0	0	0	---	---	5.208
##	12	15132141	15132141	G	T	exonic	PDE6H	nonsynonymous	SNV	PDE6H:NM_006205:exon3:c.G163T:p.G55W	12p12.3	81	10	47	0	63,8,18,2	35,0,12,0	0,0,71,10	0,0,47,0	0.123	0.8	0	---	1.877	71	10	0	47	0	0	1.877	0	0	0	0	---	---	7.481
##	12	20876048	20876048	-	G	exonic	SLCO1C1	frameshift	insertion	SLCO1C1:NM_001145944:exon7:c.692_693insG:p.L231fs,SLCO1C1:NM_001145945:exon9:c.899_900insG:p.L300fs,SLCO1C1:NM_017435:exon9:c.1046_1047insG:p.L349fs,SLCO1C1:NM_001145946:exon10:c.1046_1047insG:p.L349fs	12p12.2	97	11	57	0	82,10,15,1	50,0,7,0	---	---	0.113	0.909	0	---	2.139	76	11	0	51	0	2	2.152	0	0	1	16	---	---	6.179
```

*Not required, but you can attach "RNAseq" information.*
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

*Not required, but you can attach "Copy Number" information.*
```r
data("CopyNum")
print(CopyNum, row.names = FALSE)
```

```
##	Chromosome	Position	Log.R	segmented.LogR	BAF	segmented.BAF	Copy.number	Minor.allele	Raw.copy.number
##	1	564621	0.6071447	-0.09862298	1	NA	2	1	4.3540752
##	1	799463	0.1519967	-0.09862298	1	NA	2	1	2.1467339
##	1	1017216	0.8146911	-0.09862298	0	NA	3	1	5.8658499
##	1	1158277	-1.9594627	-0.09862298	0	NA	1	0	-0.5035897
##	1	1242215	0.2927962	-0.09862298	0	NA	2	1	2.6999875
##	1	1462766	-0.2234090	-0.09862298	1	NA	2	1	1.0726687
```

## 4. Sample Codes
lib/ccfp.jar  
lib/netMHCIIpan-3.1  
lib/netMHCpan-3.0
#
data/CopyNum.txt  
data/hla_table.txt  
data/hla_table2.txt  
data/RNAseq.txt  
data/sample.txt  
data/refFlat.cut.txt  
data/refFlat.txt  
data/refMrna.fa  
data/refMrna.merge.cut1.fa  
data/refMrna.merge.cut2.fa  
data/refMrna.merge.cut3.fa  
data/refMrna.merge.fa
##

Calculate Neoantigens on SNVs for HLA Class I. 
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
              refFlat_file = "data/refFlat.txt",
              refMrna_1 = "data/refMrna.cut1.fa",
              refMrna_3 = "data/refMrna.cut3.fa",
              RNAseq_file = "data/RNAseq.txt",
              CNV = "CopyNum.txt",
              Purity = 0.8,
              samtools_dir = "samtools-0.1.19/samtools",
              ccfp_dir = "lib/ccfp.jar",
              netMHCpan_dir = "lib/netMHCpan-3.0/netMHCpan",
              refDNA = NA)
```

Marge Results for HLA Class I. 
```
MainMergeClass1(hmdir = getwd(),
                input_dir = "data",
                input_file_prefix = "sample",
                Tumor_RNA_BASED_ON_DNA = TRUE)
```

Calculate Neoantigens for HLA Class II. 
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
              refFlat_file = "data/refFlat.txt",
              refMrna_1 = "data/refMrna.cut1.fa",
              refMrna_3 = "data/refMrna.cut3.fa",
              RNAseq_file = "data/RNAseq.txt",
              CNV = "CopyNum.txt",
              Purity = 0.8,
              samtools_dir = "samtools-0.1.19/samtools",
              ccfp_dir = "lib/ccfp.jar",
              netMHCpan_dir = "lib/netMHCIIpan-3.1/netMHCIIpan",
              refDNA = NA)
```

Marge Results for HLA Class II. 
```
MainMergeClass2(hmdir = getwd(),
                input_dir = "data",
                input_file_prefix = "sample",
                Tumor_RNA_BASED_ON_DNA = TRUE)
```

Calculate Neoantigens on Indels for HLA Class I. 
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
                refFlat_file = "data/refFlat.txt",
                refMrna_1 = "data/refMrna.cut1.fa",
                refMrna_3 = "data/refMrna.cut3.fa",
                RNAseq_file = "data/RNAseq.txt",
                CNV = NA,
                Purity = 0.8,
                samtools_dir = "samtools-0.1.19/samtools",
                ccfp_dir = "lib/ccfp.jar",
                netMHCpan_dir = "lib/netMHCpan-3.0/netMHCpan",
                refDNA = NA)
```

Calculate Neoantigens on Indels for HLA Class II. 
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
                refFlat_file = "data/refFlat.txt",
                refMrna_1 = "data/refMrna.cut1.fa",
                refMrna_3 = "data/refMrna.cut3.fa",
                RNAseq_file = "data/RNAseq.txt",
                CNV = NA,
                Purity = 0.8,
                samtools_dir = "samtools-0.1.19/samtools",
                ccfp_dir = "lib/ccfp.jar",
                netMHCpan_dir = "lib/netMHCIIpan-3.1/netMHCIIpan",
                refDNA = NA)
```

## 5. Result
Samples of result files are available at https://github.com/hase62/Neoantimon/raw/master/data/Result. 

sample.CLASS1.ALL.txt
```
##	HLA	Pos	Gene	MutatedPeptide	Mut_IC50	Mut_Rank	Norm_Peptide	Norm_IC50	Norm_Rank		Gene ID	Chr	NM_ID	Change	ref	alt	Prob	Mutation Prob.	Exon Start	Exon End	Mutation Position	Depth	TumorDepth	Peptide Normal	Peptide Mutation	TotalRNA	TumorRNARatio	TumorRNA	nA	nB	Checker	MutRatio	MutRatio Min	MutRatio Max
##	HLA-A*02:01	2	0_CYP4A11	HQERCREEIHSLP	37362.9	55.00	HQERCREEIHSLL	27284.3	32.00	1	0_CYP4A11	1	NM_000778	c.T1064C	A	G	0	0	47394845	47407156	1_47399872	49	113	KHQERCREEIHSLLGDGASITWNHLDQ	KHQERCREEIHSLPGDGASITWNHLDQ	NA	NA	NA	NA	NA	NA	NA	NA	NA
##	HLA-A*02:01	3	0_CYP4A11	QERCREEIHSLPG	45831.8	90.00	QERCREEIHSLLG	44161.6	85.00	1	0_CYP4A11	1	NM_000778	c.T1064C	A	G	0	0	47394845	47407156	1_47399872	49	113	KHQERCREEIHSLLGDGASITWNHLDQ	KHQERCREEIHSLPGDGASITWNHLDQ	NA	NA	NA	NA	NA	NA	NA	NA	NA
##	HLA-A*02:01	4	0_CYP4A11	ERCREEIHSLPGD	45989.7	90.00	ERCREEIHSLLGD	44575.4	85.00	1	0_CYP4A11	1	NM_000778	c.T1064C	A	G	0	0	47394845	47407156	1_47399872	49	113	KHQERCREEIHSLLGDGASITWNHLDQ	KHQERCREEIHSLPGDGASITWNHLDQ	NA	NA	NA	NA	NA	NA	NA	NA	NA
##	HLA-A*02:01	5	0_CYP4A11	RCREEIHSLPGDG	46482.0	95.00	RCREEIHSLLGDG	44377.6	85.00	1	0_CYP4A11	1	NM_000778	c.T1064C	A	G	0	0	47394845	47407156	1_47399872	49	113	KHQERCREEIHSLLGDGASITWNHLDQ	KHQERCREEIHSLPGDGASITWNHLDQ	NA	NA	NA	NA	NA	NA	NA	NA	NA
##	HLA-A*02:01	6	0_CYP4A11	CREEIHSLPGDGA	42936.1	80.00	CREEIHSLLGDGA	41760.9	70.00	1	0_CYP4A11	1	NM_000778	c.T1064C	A	G	0	0	47394845	47407156	1_47399872	49	113	KHQERCREEIHSLLGDGASITWNHLDQ	KHQERCREEIHSLPGDGASITWNHLDQ	NA	NA	NA	NA	NA	NA	NA	NA	NA
```

sample.CLASS1.IC50min.txt
```
##	HLA	Pos	Gene	MutatedPeptide	Mut_IC50	Mut_Rank	Norm_Peptide	Norm_IC50	Norm_Rank		Gene ID	Chr	NM_ID	Change	ref	alt	Prob	Mutation Prob.	Exon Start	Exon End	Mutation Position	Depth	TumorDepth	Peptide Normal	Peptide Mutation	TotalRNA	TumorRNARatio	TumorRNA	nA	nB	Checker	MutRatio	MutRatio Min	MutRatio Max
##	HLA-A*02:01	12	0_CYP4A11	SLPGDGASI	1217.8	4.50	SLLGDGASI	123.5	1.20	1	0_CYP4A11	1	NM_000778	c.T1064C	A	G	0	0	47394845	47407156	1_47399872	49	113	KHQERCREEIHSLLGDGASITWNHLDQ	KHQERCREEIHSLPGDGASITWNHLDQ	NA	NA	NA	NA	NA	NA	NA	NA	NA
##	HLA-A*02:01	7	1_DHX15	YLEAAIRAV	36.1	0.50	YLEAAIRTV	41.7	0.60	2	1_DHX15	4	NM_001358	c.A1012G	T	C	0	0	24529087	24586184	4_24556416	151	294	PEPERDYLEAAIRTVIQIHMCEEEEGD	PEPERDYLEAAIRAVIQIHMCEEEEGD	NA	NA	NA	NA	NA	NA	NA	NA	NA
##	HLA-A*02:01	14	2_GAS1	HCNLALSRYLTYC	14427.0	17.00	RCNLALSRYLTYC	12017.5	15.00	3	2_GAS1	9	NM_002048	c.G533A	C	T	0	0	89559276	89562104	9_89561162	26	46	GCTEARRRCDRDSRCNLALSRYLTYCG	GCTEARRRCDRDSHCNLALSRYLTYCG	NA	NA	NA	NA	NA	NA	NA	NA	NA
```





