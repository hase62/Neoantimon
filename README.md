## A manuscript is at bioarxiv (https://www.biorxiv.org/content/10.1101/869388v1).

## 0. Preliminary Use
### -Install R (Required)
CentOS
```
[command line]
yum -y install R
```

Mac(Devian)
```
[command line]
brew cask install r
```

### -Generate output without calculation
This code is a simple sample code for preliminary use to confirm the output. 

```
[R]
install.packages("devtools");
library(devtools);
install_github('hase62/Neoantimon');
library(Neoantimon);

data("sample_vcf.annovar")
data("sample_hla_table_c1")
data("sample_refFlat.grch37")
data("sample_refMrna.grch37.fa")
data("sample_result_SNV_CLASS1_ALL")

MainSNVClass1(input_annovar_format_file = sample_vcf.annovar,
              hla_types = sample_hla_table_c1[1,-1],
              refflat_file = sample_refFlat.grch37,
              refmrna_file = sample_refMrna.grch37.fa,
              netMHCpan_dir = NA)

write.table(file = "result.ID.SNV1/data.ID_SNV.peptide.SNV_CLASS1_ALL.txt", 
            x = sample_result_SNV_CLASS1_ALL[grep("0_DHX15", sample_result_SNV_CLASS1_ALL$Gene), ], 
            row.names = FALSE, quote = FALSE, sep = "\t")
```

You can get the following outputs. 

data.ID_SNV.peptide.txt (annotation file)
```
##  Number  GeneSymbol  NM_id   AAchanges   Ref Alt Prob    Mutation_Prob   Exon_Start  Exon_End    Mutation_Position   Evaluated_Mutant_Peptide    Evaluated_Wt_Peptide    ..
## 1	0_DHX15	 4	NM_001358	c.A1012G	T	C	0	0	24529097	24586177	24556416	294	143	PEPERDYLEAAIRTVIQIHMCEEEEGD	PEPERDYLEAAIRAVIQIHMCEEEEGD ..
```

data.ID_SNV.peptide.fasta (input file for netMHCpan)
```
##  >0_DHX15
##  PEPERDYLEAAIRAVIQIHMCEEEEGD
```

data.ID_SNV.wtpeptide.fasta (input file for netMHCpan)
```
##  >0_DHX15
##  PEPERDYLEAAIRTVIQIHMCEEEEGD
```

data.ID_SNV.peptide.SNV_CLASS1_ALL.txt (output file)
```
##  HLA	Pos	Gene	Evaluated_Mutant_Peptide	Mut_IC50	Mut_Rank	Evaluated_Wt_Peptide	Wt_IC50	Wt_Rank	Chr	NM_ID	Change	Ref	Alt	Prob	Mutation_Prob.	Exon_Start	Exon_End	Mutation_Position	Total_Depth	Tumor_Depth	Wt_Peptide	Mutant_Peptide	Total_RNA	Tumor_RNA_Ratio	Tumor_RNA	Tumor_RNA_based_on_DNA	MutRatio	MutRatio_Min	MutRatio_Max
##  HLA-A*02:01	2	0_DHX15	EPERDYLEAAIRA	37952.9	63.3899	EPERDYLEAAIRT	41213.3	75.7538	4	NM_001358	c.A1012G	T	C	0	0	24529087	24586184	4_24556416	294	143	PEPERDYLEAAIRTVIQIHMCEEEEGD	PEPERDYLEAAIRAVIQIHMCEEEEGD	1.35204	NA	NA	0.657624897959184	NA	NA	NA
##  HLA-A*02:01	3	0_DHX15	PERDYLEAAIRAV	12108.2	16.51	PERDYLEAAIRTV	11859.7	16.2523	4	NM_001358	c.A1012G	T	C	0	0	24529087	24586184	4_24556416	294	143	PEPERDYLEAAIRTVIQIHMCEEEEGD	PEPERDYLEAAIRAVIQIHMCEEEEGD	1.35204	NA	NA	0.657624897959184	NA	NA	NA
##  HLA-A*02:01	4	0_DHX15	ERDYLEAAIRAVI	12109.5	16.5114	ERDYLEAAIRTVI	14831.9	19.2905	4	NM_001358	c.A1012G	T	C	0	0	24529087	24586184	4_24556416	294	143	PEPERDYLEAAIRTVIQIHMCEEEEGD	PEPERDYLEAAIRAVIQIHMCEEEEGD	1.35204	NA	NA	0.657624897959184	NA	NA	NA
##  ..
##  
```

## 1. Preparation
### -Install wget (Required)
CentOS
```
[command line]
yum install wget
```

Mac(Devian)
```
[command line]
brew install wget
```

### -Download and Set netMHCpan4.0 (Required)

1. Download netMHCpan4.0 from https://services.healthtech.dtu.dk/service.php?NetMHCpan-4.0 and move it to any working directory. 
(We assume that "lib" directory contains netMHCpan-4.0a.{Darwin|Linux}.tar.gz.)

2. Run the initial setting script at the directory that contains downloaded  as followings.

We assume
```
lib/
  ├ netMHCpan-4.0a.{Darwin|Linux}.tar
  └ setNetMHCpan4.0.sh
```

Run
```
[command line]
wget --no-check-certificate https://github.com/hase62/Neoantimon/raw/master/lib/setNetMHCpan4.0.sh
chmod 750 setNetMHCpan4.0.sh
./setNetMHCpan4.0.sh 
```

We have
```
lib/
    ├ netMHCpan-4.0a.{Darwin|Linux}.tar
    ├ setNetMHCpan4.0.sh
    └ NetMHCpan4.0/
      ├ {Darwin|Linux}_x86_64
      ├ data
      ├ data.{Darwin|Linux}.tar.gz
      ├ netMHCpan
      ├ netMHCpan-4.0.readme
      ├ netMHCpan-e
      ├ netMHCpan.1
      ├ test
      └ tmp
```

### -Download and Set mhcflurry (Not Required)
1. (Recommended) Install anaconda from https://www.anaconda.com/distribution/, and then run the following codes. 
```
[command line]
pip install mhcflurry
mhcflurry-downloads fetch
pip install mhctools
```

2. Otherwise, install python from https://www.python.org/downloads/release, and then run the above codes.

### -Download and Set netMHCIIpan3.2 (Required)

1. Download netMHCIIpan 3.2 from https://services.healthtech.dtu.dk/service.php?NetMHCIIpan-3.2 and move it to the working directory. 
(We assume that "lib" directory contains netMHCIIpan-3.2.{Darwin|Linux}.tar.gz.)

2. Do initial setting at the working directory as followings.

We assume
```
lib/
  ├ netMHCIIpan-3.2.{Darwin|Linux}.tar
  └ setNetMHCIIpan3.2.sh
```

Run
```
[command line]
wget --no-check-certificate https://github.com/hase62/Neoantimon/raw/master/lib/setNetMHCIIpan3.2.sh
chmod 750 setNetMHCIIpan3.2.sh
./setNetMHCIIpan3.2.sh
```

We have
```
lib/
    ├ netMHCIIpan-3.2.{Darwin|Linux}.tar
    ├ setNetMHCIIpan3.2.sh
    └ netMHCIIpan-3.2/
      ├ {Darwin|Linux}_x86_64
      ├ data
      ├ data.{Darwin|Linux}.tar.gz
      ├ netMHCIIpan
      ├ NetMHCIIpan-3.2.pl
      ├ netMHCIIpan-3.2.readme
      ├ netMHCIIpan-e
      ├ etMHCIIpan.1
      ├ test
      └ tmp
```

### -Download refMrna Files (Required)
(You have to select your corresponding version from GRCh38, hg38, GRCh37 or hg19.)

**GRCh38/hg38**: Run the following codes. 
```
[command line]
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/refMrna.fa.gz
gunzip refMrna.fa.gz
mv refMrna.fa refMrna.grch38.fa
```

**GRCh37/hg19**: Run the following codes. 
```
[command line]
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/refMrna.fa.gz
gunzip refMrna.fa.gz
mv refMrna.fa refMrna.grch37.fa
```

### -Download refFlat Files (Required)
(You have to select your corresponding version from GRCh38, hg38, GRCh37 or hg19.)

**GRCh38/hg38**: Run the following codes. 
```
[command line]
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/refFlat.txt.gz
gunzip refFlat.txt.gz
mv refFlat.txt refFlat.grch38.txt
```

**GRCh37/hg19**: Run the following codes. 
```
[command line]
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/refFlat.txt.gz
gunzip refFlat.txt.gz
mv refFlat.txt refFlat.grch37.txt
```

### -Install Samtools (Required to analyze structural variants or utilize RNA BAM)

1. (Recommended) Install anaconda from https://www.anaconda.com/distribution/, and then run the following codes. 
```
[command line]
conda install -c bioconda samtools
conda install -c bioconda/label/cf201911 samtools
```

2. Otherwise, you can install local samtools as followings. 
```
[command line]
wget https://sourceforge.net/projects/samtools/files/samtools/0.1.19/samtools-0.1.19.tar.bz2
tar jxf samtools-0.1.19.tar.bz2
```

### -Download human refSeq (Required to analyze structural variants or utilize RNA BAM)

(You have to select your corresponding version from GRCh38, hg38, GRCh37 or hg19.)

**GRCh38**: Run the following codes.
```
[command line]
wget ftp://ftp.ensembl.org/pub/release-87/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.toplevel.fa.gz
mv Homo_sapiens.GRCh38.dna.toplevel.fa.gz GRCh38.fa.gz
gunzip GRCh38.fa.gz
samtools faidx GRCh38.fa
```

**hg38**: Run the following codes.
```
[command line]
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
gunzip hg38.fa.gz
samtools faidx hg38.fa
```

**GRCh37/hg19**: Run the following codes.
```
[command line]
wget ftp://ftp.ensembl.org/pub/release-75//fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.75.dna.toplevel.fa.gz
mv Homo_sapiens.GRCh37.75.dna.toplevel.fa.gz GRCh37.fa.gz
gunzip GRCh37.fa.gz
samtools faidx GRCh37.fa
```

## 2. Prepare to Use on R
**Required**
```
[R]
install.packages("devtools");
library(devtools);
install_github('hase62/Neoantimon');
library(Neoantimon);
```

**Suggested to Install for Reading Data at High Speed**
```
[R]
install.packages('data.table');
library(data.table);
```

## 3. Data Format
### -HLA Table (or, You can indicate HLA types directly)

**1. A HLA Class I table file must be according to the following format.**

```r
[R]
data("sample_hla_table_c1")
print(sample_hla_table_c1, row.names = FALSE)
```

```
##	Name	A1	A2	B1	B2	C1	C2
##	sample	A*02:01	A*32:01	B*15:17	B*51:01	C*07:01	C*15:02
##	sample2	A*02:01	A*32:01	B*15:17	B*51:01	C*07:01	C*15:02
```

**2. A HLA Class II table file must be according to the following format.**

```r
[R]
data("sample_hla_table_c2")
print(sample_hla_table_c2, row.names = FALSE)
```

```
##	Name	DPA11	DPA12	DPB11	DPB12	DQA11	DQA12	DQB11	DQB12	DRB11	DRB12
##	sample	DPA1*01:03	DPA1*02:01	DPB1*02:01	DPB1*09:01	DQA1*01:02	DQA1*05:05	DQB1*03:01	DQB1*06:04	DRB1*11:04	DRB1*13:02
##	sample2	DPA1*01:03	DPA1*02:01	DPB1*02:01	DPB1*09:01	DQA1*01:02	DQA1*05:05	DQB1*03:01	DQB1*06:04	DRB1*11:04	DRB1*13:02
```

### -Annotated vcf file
**An annovar format or an Ensembl Variant Effect Predictor (VEP) format annotated vcf file is required for Snv/Indel.**

#### -Annotated vcf file by Annovar
When indicating annovar format, it must include columns representing "Chromosome Number", "Mutation Start Position", "Mutation End Position", "Mutation Ref", "Mutation Alt", and "NM_ID (AAChange.refGene)".
Annotations "Chr", "Start", "End", "Ref", "Alt", "AAChange.refGene", "Depth_tumor", and "Depth_normal" are automatically detected. Otherwise, you have to manually indicate columns. 
```r
[R]
data("sample_vcf.annovar")
print(sample_vcf.annovar, row.names = FALSE)
```

```
## Chr     Start       End Ref Alt Func.refGene Gene.refGene   GeneDetail.refGene ExonicFunc.refGene                              AAChange.refGene   cytoBand depth_tumor variantNum_tumor depth_normal
##   1 116941338 116941338   T   C       exonic       ATP1A1           synonymous                SNV   ATP1A1:NM_001160234:exon16:c.T2127C:p.D709D     1p13.1         100               39          111
##   4  24556416  24556416   T   C       exonic        DHX15        nonsynonymous                SNV        DHX15:NM_001358:exon5:c.A1012G:p.T338A     4p15.2         143               47          151
##   4  70156404  70156404   -   T       exonic      UGT2B28           frameshift          insertion   UGT2B28:NM_053039:exon5:c.1186dupT:p.L395fs     4q13.2          43               15           41
##   6  75899298  75899298   T   -       exonic      COL12A1           frameshift           deletion    COL12A1:NM_004370:exon6:c.628delA:p.I210fs       6q13         122               38           73
##   9  89561162  89561162   C   T       exonic         GAS1        nonsynonymous                SNV          GAS1:NM_002048:exon1:c.G533A:p.R178H    9q21.33          20                5           26
...
```

#### -Annotated file by VEP
```r
[R]
data("sample_vcf.vep")
print(sample_vcf.vep, row.names = FALSE)
```

```
##  Uploaded_variation	Location	Allele	Gene	Feature	Feature_type	Consequence	cDNA_position	CDS_position	Protein_position	Amino_acids	Codons	Existing_variation	Extra
##  file_name	1:1290161	A	ENSG00000162576	ENST00000309212	Transcript	missense_variant	881	850	284	R/C	Cgc/Tgc	-	IMPACT=MODERATE;STRAND=-1
##  file_name	1:1434364	T	ENSG00000160072	ENST00000474481	Transcript	downstream_gene_variant	-	-	-	-	-	-	IMPACT=MODIFIER;DISTANCE=2783;STRAND=1
##  file_name	1:12919618	A	ENSG00000120952	ENST00000240189	Transcript	missense_variant	445	358	120	W/R	Tgg/Agg	-	IMPACT=MODERATE;STRAND=1
##  file_name	1:16890602	G	ENSG00000219481	ENST00000430580	Transcript	missense_variant	4144	3256	1086	M/L	Atg/Ctg	-	IMPACT=MODERATE;STRAND=-1
##  file_name	1:16918457	C	ENSG00000219481	ENST00000392963	Transcript	"missense_variant,NMD_transcript_variant"	518	60	20	I/M	atC/atG	-	IMPACT=MODERATE;STRAND=-1
##  file_name	1:20297876	G	ENSG00000188257	ENST00000482011	Transcript	downstream_gene_variant	-	-	-	-	-	-	IMPACT=MODIFIER;DISTANCE=4290;STRAND=-1
##  file_name	1:48906934	T	ENSG00000132122	ENST00000487543	Transcript	"intron_variant,NMD_transcript_variant"	-	-	-	-	-	-	IMPACT=MODIFIER;STRAND=-1
...
```

Please install the following library when you use vep-annotated files. 
```
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("biomaRt")
library(biomaRt)
```

For example, please run as following. 
```
  MainSNVClass1(input_vep_format_file = sample_vcf.vep,
                hla_types = sample_hla_table_c1[1,-1],
                refflat_file = sample_refFlat.grch37,
                refmrna_file = sample_refMrna.grch37.fa,
                netMHCpan_dir = NA)
```

### -Annotated BND format vcf file

*An annotated BND format vcf file is required for SV fusion.*

It must include columns representing "Chromosome Number", "Mutation Start Position", "Mutation End Position", "Mutation Ref", "Mutation Alt", and "NM_ID (AAChange.refGene)" or "Gene Symbol (Gene.refGene)".
Annotations "Chr", "Start", "End", "Ref", "Alt", "Depth_tumor", and "Depth_normal" are automatically detected. Otherwise, you have to manually indicate columns. 
```r
[R]
data("sample_sv_bnd")
print(sample_sv_bnd, row.names = FALSE)
```

```
## Chr     Start       End Ref            Alt Func.refGene Gene.refGene           mateID
##   1 115005805 115005805   C C]20:34827929]       exonic       TRIM33     SVMERGE137_1
##   1 204908711 204908711   A [2:172743385[A     intronic        NFASC      SVMERGE15_1
##   2  25870534  25870534   T  ]2:25965720]G     intronic         DTNB     SVMERGE116_1
##   2  25965720  25965720   T  T[2:25870536[       exonic        ASXL2     SVMERGE116_2
##   2 214794791 214794791   C C[2:214798169[       exonic       SPAG16       SVMERGE3_1
##   2 214798169 214798169   T ]2:214794791]T     intronic       SPAG16       SVMERGE3_2
...
```

### -RNA expression file

An RNA expressoin file is not required, but you can attach "RNA expression" information by indicating "rnaexp_file".
If you also indicate "rnabam_file", variant allele frequencies and tumor specific RNA expressions are also attached to the results. 

```r
[R]
data("sample_rna_exp")
print(sample_rna_exp, row.names = FALSE)
```

```
##  gene_short_name		   ChromosomeNum:Start-Endlocus expression
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

### -CNV file

A copynumber file is not required, but you can attach "Copy Number" information by indicating "cnv_file" and "purity" in main functions.
They are used to calculate tumor subclonal cell population. 
Purity is set 1 as default value. 

```r
[R]
data("sample_copynum")
print(sample_copynum, row.names = FALSE)
```

```
##	Chromosome	Position	     Log.R	segmented.LogR	BAF	segmented.BAF	Copy.number	Minor.allele	Raw.copy.number
##			 1	  564621	 0.6071447	 -0.09862298	  1			   NA			  2			   1		  4.3540752
##			 1	  799463	 0.1519967	 -0.09862298	  1			   NA			  2			   1		  2.1467339
##			 1	  1017216	 0.8146911	 -0.09862298	  0			   NA			  3			   1		  5.8658499
##			 1	  1158277	-1.9594627	 -0.09862298	  0			   NA			  1			   0		 -0.5035897
##			 1	  1242215	 0.2927962	 -0.09862298	  0			   NA			  2			   1		  2.6999875
##			 1	  1462766	-0.2234090	 -0.09862298	  1			   NA			  2			   1		  1.0726687
```

## 4. Sample Codes

### -Download SampleFiles

Run the following codes. 
(We assume that data will be downloaded into "lib" directory)
```
[command line]
wget --no-check-certificate https://github.com/hase62/Neoantimon/raw/master/lib/data.zip
unzip data.zip
```

We assume the following directory structure according to 1. Preparation.
```
~/opt/anaconda3
lib/
    ├ NetMHCpan4.0
    ├ netMHCIIpan-3.2
    ├ refFlat.grch37.txt
    ├ refMrna.grch37.fa
    └ data/
      ├ sample_result_INDEL_CLASS1_ALL.txt
      ├ sample_result_INDEL_CLASS2_ALL.txt
      ├ sample_result_SeqFragment_CLASS1_ALL.txt
      ├ sample_result_SeqFragment_CLASS2_ALL.txt
      ├ sample_result_SNV_CLASS1_ALL.txt
      ├ sample_result_SNV_CLASS2_ALL.txt
      ├ sample_result_SVFusion_CLASS1_ALL.txt
      ├ sample_result_SVFusion_CLASS2_ALL.txt
      ├ sample_copynum.txt
      ├ sample_hla_table_c1.txt
      ├ sample_hla_table_c2.txt
      ├ sample_refMrna.grch37.fa.txt
      ├ sample_refFlat.grch37.txt
      ├ sample_rna_exp.txt
      ├ sample_vcf.annovar.txt
      ├ sample_vcf.vep.txt
      ├ sample_sv_bnd.txt
      └ sample.snps.vcf
```

#### Calculate Neoantigens on SNVs for HLA Class I and II. 
<kbd><img src="https://github.com/hase62/Neoantimon/blob/images/images/ForExplanation_snv.png" width="640px"></kbd>

To calculate the binding affinity of neoantigen candaites, which are generated from from SNVs, to HLA ClassI. 
When using MHCflurry, [[1]] and [[2]] include the results of NetMHCpan and MHCflurry, respectively. 
```
[R]
  Result_HLA1_SNV <- MainSNVClass1(input_annovar_format_file = "data/sample_vcf.annovar.txt",
                                   file_name_in_hla_table = "sample",
                                   hla_file = "data/sample_hla_table_c1.txt",
                                   refflat_file  = "refFlat.grch37.txt",
                                   refmrna_file = "refMrna.grch37.fa",
                                   rnaexp_file = "data/sample_rna_exp.txt",
                                   netMHCpan_dir = "netMHCpan-4.0/netMHCpan",
                                   depth_tumor_column = 12,
                                   depth_normal_column = 14,
                                   SNPs = "data/sample.snps.vcf",
                                   multiple_variants = TRUE,
                                   MHCflurry = "~/opt/anaconda3/bin/mhctools")
```

To calculate the binding affinity of neoantigen candaites, which are generated from from SNVs, to HLA ClassII. 
```
[R]
  Result_HLA2_SNV <- MainSNVClass2(input_annovar_format_file = "data/sample_vcf.annovar.txt",
                                   file_name_in_hla_table = "sample",
                                   hla_file = "data/sample_hla_table_c2.txt",
                                   refflat_file  = "refFlat.grch37.txt",
                                   refmrna_file = "refMrna.grch37.fa",
                                   rnaexp_file = "data/sample_rna_exp.txt",
                                   netMHCIIpan_dir = "netMHCIIpan-3.2/netMHCIIpan",
                                   depth_tumor_column = 12,
                                   depth_normal_column = 14,
                                   SNPs = "data/sample.snps.vcf",
                                   multiple_variants = TRUE)
```

#### Calculate Neoantigens on INDELs for HLA Class I and II. 
<kbd><img src="https://github.com/hase62/Neoantimon/blob/images/images/ForExplanation_indel.png" width="640px"></kbd>

To calculate the binding affinity of neoantigen candaites, which are generated from from indels, to HLA ClassI. 
When using MHCflurry, [[1]] and [[2]] include the results of NetMHCpan and MHCflurry, respectively. 
```
[R]
  Result_HLA1_INDEL <- MainINDELClass1(input_annovar_format_file = "data/sample_vcf.annovar.txt",
                                       file_name_in_hla_table = "sample",
                                       hla_file = "data/sample_hla_table_c1.txt",
                                       refflat_file  = "refFlat.grch37.txt",
                                       refmrna_file = "refMrna.grch37.fa",
                                       rnaexp_file = "data/sample_rna_exp.txt",
                                       netMHCpan_dir = "netMHCpan-4.0/netMHCpan",
                                       depth_tumor_column = 12,
                                       depth_normal_column = 14,
                                       SNPs = "data/sample.snps.vcf",
                                       multiple_variants = TRUE,
                                       MHCflurry = "~/opt/anaconda3/bin/mhctools")
```

To calculate the binding affinity of neoantigen candaites, which are generated from from indels, to HLA ClassII. 
```
[R]
  Result_HLA2_INDEL <- MainINDELClass2(input_annovar_format_file = "data/sample_vcf.annovar.txt",
                                       file_name_in_hla_table = "sample",
                                       hla_file = "data/sample_hla_table_c2.txt",
                                       refflat_file  = "refFlat.grch37.txt",
                                       refmrna_file = "refMrna.grch37.fa",
                                       rnaexp_file = "data/sample_rna_exp.txt",
                                       netMHCIIpan_dir = "netMHCIIpan-3.2/netMHCIIpan",
                                       depth_tumor_column = 12,
                                       depth_normal_column = 14,
                                       SNPs = "data/sample.snps.vcf",
                                       multiple_variants = TRUE)
```

#### Calculate Neoantigens on SV fusions for HLA Class I and II. 
<kbd><img src="https://github.com/hase62/Neoantimon/blob/images/images/ForExplanation_sv.png" width="640px"><kbd>

To calculate the binding affinity of neoantigen candaites, which are generated from from SVs, to HLA ClassI. 
```
[R]
  Result_HLA1_SV <- MainSVFUSIONClass1(input_file = "data/sample_sv_bnd.txt",
                                       file_name_in_hla_table = "sample",
                                       hla_file = "data/sample_hla_table_c1.txt",
                                       refflat_file  = "refFlat.grch37.txt",
                                       refmrna_file = "refMrna.grch37.fa",
                                       rnaexp_file = "data/sample_rna_exp.txt",
                                       netMHCpan_dir = "netMHCpan-4.0/netMHCpan",
                                       refdna_file = "GRCh37.fa",
                                       mutation_alt_bnd_column = 5,
                                       gene_symbol_column = 7,
                                       mate_id_column = 8)
```

To calculate the binding affinity of neoantigen candaites, which are generated from from SVs, to HLA ClassII. 
```
[R]
  Result_HLA2_SV <- MainSVFUSIONClass2(input_file = "data/sample_sv_bnd.txt",
                                       file_name_in_hla_table = "sample",
                                       hla_file = "data/sample_hla_table_c2.txt",
                                       refflat_file  = "refFlat.grch37.txt",
                                       refmrna_file = "refMrna.grch37.fa",
                                       rnaexp_file = "data/sample_rna_exp.txt",
                                       netMHCIIpan_dir = "netMHCIIpan-3.2/netMHCIIpan",
                                       refdna_file = "GRCh37.fa",
                                       mutation_alt_bnd_column = 5,
                                       gene_symbol_column = 7,
                                       mate_id_column = 8)
```

#### Calculate Neoantigens from a fragment of RNA sequence for HLA Class I and II by comparing to the original protein. 
<kbd><img src="https://github.com/hase62/Neoantimon/blob/images/images/ForExplanation_rna.png" width="640px"></kbd>

To calculate the binding affinity of neoantigen candaites, which are directly generated from RNA sequences, to HLA ClassI. 
The peptides included in the original genes ("NM_003998", "NM_001165412") are removed from the results. 
```
[R]
  Result_HLA1_Seq <- MainSeqFragmentClass1(input_sequence = "atggcagaagatgatccatatttgggaaggcctgaaaaaatgtttcatttggatccttctttgactcatacaatatttaatc",
                                           file_name_in_hla_table = "sample",
                                           hla_file = "data/sample_hla_table_c1.txt",
                                           hmdir = getwd(),
                                           job_id = "NO_job_id",
                                           refflat_file  = "refFlat.grch37.txt",
                                           refmrna_file = "refMrna.grch37.fa",
                                           netMHCpan_dir = "netMHCpan-4.0/netMHCpan",
                                           reference_nm_id = c("NM_003998", "NM_001165412"))
```

To calculate the binding affinity of neoantigen candaites, which are directly generated from RNA sequences, to HLA ClassII. 
The peptides included in the riginal genes ("NFKB1", "BCL3") are removed from the results. 
```
[R]
  Result_HLA2_Seq <- MainSeqFragmentClass2(input_sequence = "atggcagaagatgatccatatttgggaaggcctgaacaaatgtttcatttgatccttctttgactcatacaatatttaatc",
                                           file_name_in_hla_table = "sample",
                                           hla_file = "data/sample_hla_table_c2.txt",
                                           hmdir = getwd(),
                                           job_id = "NO_job_id",
                                           refflat_file  = "refFlat.grch37.txt",
                                           refmrna_file = "refMrna.grch37.fa",
                                           netMHCIIpan_dir = "netMHCIIpan-3.2/netMHCIIpan",
                                           reference_gene_symbol = c("NFKB1", "BCL3"))
```

#### Attach Opitional Information
Users can optionally provide (a) RNA expression data with and without (b) the corresponding Binary Alignment Map (BAM) file to attach total and allele specific RNA expression levels, and (c) copy number variation data to calculate the posterior probability distribution over cancer-cell fraction (CCFP) to evaluate tumor sub-clonality (Lohr et al., 2014).

<kbd><img src="https://github.com/hase62/Neoantimon/blob/images/images/ForExplanation_op2.png" width="640px"></kbd>

In addition, users can consider specific cases of existing SNPs on the mutant peptide by providing (c) SNPs data, and multiple SNVs on the same mutant peptides and among the frameshift region caused by indels. These cases are explained as followings.

<kbd><img src="https://github.com/hase62/Neoantimon/blob/images/images/ForExplanation_op1.png" width="640px"></kbd>

## 5. Result
#### Output Result

Results generated from SNVs. 
```
[R]
print(head(Result_HLA1_SNV[[1]]))
```
  
```
##	HLA	Pos	Gene	Evaluated_Mutant_Peptide	Mut_IC50	Mut_Rank	Evaluated_Wt_Peptide	Wt_IC50	Wt_Rank	Chr	NM_ID	Change	Ref	Alt	Prob	Mutation_Prob.	Exon_Start	Exon_End	Mutation_Position	Total_Depth	Tumor_Depth	Wt_Peptide	Mutant_Peptide	Total_RNA	Tumor_RNA_Ratio	Tumor_RNA	Tumor_RNA_based_on_DNA	MutRatio	MutRatio_Min	MutRatio_Max
##	HLA-A*02:01	2	0_CYP4A11	HQERCREEIHSLP	38277.9	64.5038	HQERCREEIHSLL	21979.5	27.7445	1	NM_000778	c.T1064C	A	G	0	0	47394859	47407148	1_47399872	113	64	KHQERCREEIHSLLGDGASITWNHLDQ	KHQERCREEIHSLPGDGASITWNHLDQ	0	NA	NA	0	NA	NA	NA
##	HLA-A*02:01	3	0_CYP4A11	QERCREEIHSLPG	44731.5	90.084	QERCREEIHSLLG	43057.5	83.3295	1	NM_000778	c.T1064C	A	G	0	0	47394859	47407148	1_47399872	113	64	KHQERCREEIHSLLGDGASITWNHLDQ	KHQERCREEIHSLPGDGASITWNHLDQ	0	NA	NA	0	NA	NA	NA
##	HLA-A*02:01	4	0_CYP4A11	ERCREEIHSLPGD	45387.7	92.3912	ERCREEIHSLLGD	43941.4	86.9293	1	NM_000778	c.T1064C	A	G	0	0	47394859	47407148	1_47399872	113	64	KHQERCREEIHSLLGDGASITWNHLDQ	KHQERCREEIHSLPGDGASITWNHLDQ	0	NA	NA	0	NA	NA	NA
##	HLA-A*02:01	5	0_CYP4A11	RCREEIHSLPGDG	45916.2	94.2252	RCREEIHSLLGDG	43591.4	85.5092	1	NM_000778	c.T1064C	A	G	0	0	47394859	47407148	1_47399872	113	64	KHQERCREEIHSLLGDGASITWNHLDQ	KHQERCREEIHSLPGDGASITWNHLDQ	0	NA	NA	0	NA	NA	NA
##	HLA-A*02:01	6	0_CYP4A11	CREEIHSLPGDGA	41190.5	75.6598	CREEIHSLLGDGA	40373.3	72.4303	1	NM_000778	c.T1064C	A	G	0	0	47394859	47407148	1_47399872	113	64	KHQERCREEIHSLLGDGASITWNHLDQ	KHQERCREEIHSLPGDGASITWNHLDQ	0	NA	NA	0	NA	NA	NA
##	HLA-A*02:01	7	0_CYP4A11	REEIHSLPGDGAS	42224.3	79.882	REEIHSLLGDGAS	41200.3	75.7003	1	NM_000778	c.T1064C	A	G	0	0	47394859	47407148	1_47399872	113	64	KHQERCREEIHSLLGDGASITWNHLDQ	KHQERCREEIHSLPGDGASITWNHLDQ	0	NA	NA	0	NA	NA	NA
```

Count the number of neoantigens of which thresholds are IC 50 of mutatnt peptides < 500 and IC 50 of wild-type peptides > 500. 
```
[R]
print(Export_Summary_SNV(Input = Result_HLA1_SNV[[1]], Mut_IC50_th = 500, Wt_IC50_th = 500))
```

```
##	Num_All_Alteration	Num_Evaluated_Alteration Num_Alteration_Generating_NeoAg	Num_All_Peptide	Num_Evaluated_Peptide	Num_Peptide_Generating_NeoAg
##	4	4	2	252	252	3
```

Results generated from SNVs. 
```
[R]
print(head(Result_HLA2_SNV))
```

```
##	HLA	Pos	Gene	Evaluated_Mutant_Peptide	Mut_IC50	Mut_Rank	Evaluated_Wt_Peptide	Wt_IC50	Wt_Rank	Chr	NM_ID	Change	Ref	Alt	Prob	Mutation_Prob.	Exon_Start	Exon_End	Mutation_Position	Total_Depth	Tumor_Depth	Wt_Peptide	Mutant_Peptide	Total_RNA	Tumor_RNA_Ratio	Tumor_RNA	Tumor_RNA_based_on_DNA	MutRatio	MutRatio_Min	MutRatio_Max
##	HLA-DPA10103-DPB10201	2	0_CYP4A11	PKHQERCREEIHSLP	12151.54	95	PKHQERCREEIHSLL	5209.84	75	1	NM_000778	c.T1064C	A	G	0	0	47394859	47407148	1_47399872	113	64	HPKHQERCREEIHSLLGDGASITWNHLDQMP	HPKHQERCREEIHSLPGDGASITWNHLDQMP	0	NA	NA	0	NA	NA	NA
##	HLA-DPA10103-DPB10201	3	0_CYP4A11	KHQERCREEIHSLPG	11550.05	90	KHQERCREEIHSLLG	4635.67	70	1	NM_000778	c.T1064C	A	G	0	0	47394859	47407148	1_47399872	113	64	HPKHQERCREEIHSLLGDGASITWNHLDQMP	HPKHQERCREEIHSLPGDGASITWNHLDQMP	0	NA	NA	0	NA	NA	NA
##	HLA-DPA10103-DPB10201	4	0_CYP4A11	HQERCREEIHSLPGD	12677.06	95	HQERCREEIHSLLGD	4880.31	75	1	NM_000778	c.T1064C	A	G	0	0	47394859	47407148	1_47399872	113	64	HPKHQERCREEIHSLLGDGASITWNHLDQMP	HPKHQERCREEIHSLPGDGASITWNHLDQMP	0	NA	NA	0	NA	NA	NA
##	HLA-DPA10103-DPB10201	5	0_CYP4A11	QERCREEIHSLPGDG	13675.38	95	QERCREEIHSLLGDG	5435.5	75	1	NM_000778	c.T1064C	A	G	0	0	47394859	47407148	1_47399872	113	64	HPKHQERCREEIHSLLGDGASITWNHLDQMP	HPKHQERCREEIHSLPGDGASITWNHLDQMP	0	NA	NA	0	NA	NA	NA
##	HLA-DPA10103-DPB10201	6	0_CYP4A11	ERCREEIHSLPGDGA	13240.19	95	ERCREEIHSLLGDGA	5173.27	75	1	NM_000778	c.T1064C	A	G	0	0	47394859	47407148	1_47399872	113	64	HPKHQERCREEIHSLLGDGASITWNHLDQMP	HPKHQERCREEIHSLPGDGASITWNHLDQMP	0	NA	NA	0	NA	NA	NA
##	HLA-DPA10103-DPB10201	7	0_CYP4A11	RCREEIHSLPGDGAS	13790.24	95	RCREEIHSLLGDGAS	5387.49	75	1	NM_000778	c.T1064C	A	G	0	0	47394859	47407148	1_47399872	113	64	HPKHQERCREEIHSLLGDGASITWNHLDQMP	HPKHQERCREEIHSLPGDGASITWNHLDQMP	0	NA	NA	0	NA	NA	NA
```

Count the number of neoantigens of which thresholds are IC 50 of mutatnt peptides < 500 and IC 50 of wild-type peptides > 500. 
```
[R]
print(Export_Summary_SNV(Input = Result_HLA2_SNV, Mut_IC50_th = 500, Wt_IC50_th = 500))
```

``` 
##	Num_All_Alteration	Num_Evaluated_Alteration	Num_Alteration_Generating_NeoAg	Num_All_Peptide	Num_Evaluated_Peptide	Num_Peptide_Generating_NeoAg	
##	4	4	2	60	60	10	
```

Results generated from indels. 
```
[R]
print(head(Result_HLA1_INDEL[[1]]))
```

```
HLA	Pos	Gene	Evaluated_Mutant_Peptide_Core	Evaluated_Mutant_Peptide	Mut_IC50	Mut_Rank	Chr	NM_ID	Change	Ref	Alt	Prob	Mutation_Prob.	Exon_Start	Exon_End	Mutation_Position	Total_Depth	Tumor_Depth	Wt_Peptide	Mutant_Peptide	Total_RNA	Tumor_RNA_Ratio	Tumor_RNA	Tumor_RNA_based_on_DNA	MutRatio	MutRatio_Min	MutRatio_Max
HLA-A*02:01	1	0_UGT2B28	GIPMVGIPLV	IYHGIPMVGIPLV	8158.7	12.6683	4	NM_053039	Out_c.1186dupT	-	T	0	0	70146192	70160768	4_70156404	84	43	IYHGIPMVGIPLFWDQPDNIAHMKAKGA	IYHGIPMVGIPLVLGSTX	0	NA	NA	0	NA	NA	NA
HLA-A*02:01	2	0_UGT2B28	YHGIPMVGIPLVL	YHGIPMVGIPLVL	264.7	1.9668	4	NM_053039	Out_c.1186dupT	-	T	0	0	70146192	70160768	4_70156404	84	43	IYHGIPMVGIPLFWDQPDNIAHMKAKGA	IYHGIPMVGIPLVLGSTX	0	NA	NA	0	NA	NA	NA
HLA-A*02:01	3	0_UGT2B28	HGIPMVGIPLVL	HGIPMVGIPLVLG	14185.4	18.6104	4	NM_053039	Out_c.1186dupT	-	T	0	0	70146192	70160768	4_70156404	84	43	IYHGIPMVGIPLFWDQPDNIAHMKAKGA	IYHGIPMVGIPLVLGSTX	0	NA	NA	0	NA	NA	NA
HLA-A*02:01	4	0_UGT2B28	GIPMVGIPLV	GIPMVGIPLVLGS	13299.6	17.694	4	NM_053039	Out_c.1186dupT	-	T	0	0	70146192	70160768	4_70156404	84	43	IYHGIPMVGIPLFWDQPDNIAHMKAKGA	IYHGIPMVGIPLVLGSTX	0	NA	NA	0	NA	NA	NA
HLA-A*02:01	5	0_UGT2B28	IPMVGIPLVL	IPMVGIPLVLGST	17674.7	22.4056	4	NM_053039	Out_c.1186dupT	-	T	0	0	70146192	70160768	4_70156404	84	43	IYHGIPMVGIPLFWDQPDNIAHMKAKGA	IYHGIPMVGIPLVLGSTX	0	NA	NA	0	NA	NA	NA
HLA-A*02:01	2	0_UGT2B28	GIPMVGIPLV	YHGIPMVGIPLV	5439.5	9.9316	4	NM_053039	Out_c.1186dupT	-	T	0	0	70146192	70160768	4_70156404	84	43	IYHGIPMVGIPLFWDQPDNIAHMKAKGA	IYHGIPMVGIPLVLGSTX	0	NA	NA	0	NA	NA	NA
```

Count the number of neoantigens of which threshold is IC 50 of mutatnt peptides < 500. 
```
[R]
print(Export_Summary_IndelSV(Input = Result_HLA1_INDEL[[1]], Mut_IC50_th = 500))
```

```
##	Num_All_Alteration	Num_Evaluated_Alteration	Num_Alteration_Generating_NeoAg	Num_All_Peptide	Num_Evaluated_Peptide	Num_Peptide_Generating_NeoAg		
##	3	3	3	270	270	15		
```

Count the number of neoantigens for each peptide of which threshold is IC 50 of mutatnt peptides < 500. 
```
[R]
print(Export_Summary_IndelSV_perFragments(Input = Result_HLA1_INDEL[[1]], Mut_IC50_th = 500))
```

```
##		EAIYHGIPMVGIPLVLGSTX-0_UGT2B28	QYYQRDELLAAIKKFHIKVATQX-1_COL12A1	IMEMARDFLPSLKNPFWKPSILPIFMYKHCSVQFSVRHGDVQTKVHX-2_SLCO1C1
##	Num_Peptide_Per_Pep	5	8	32
##	Num_Cond_Peptide_Per_Pep	5	8	32
##	Num_Rest_Peptide_Per_Pep	5	7	20
##	Num_Rest_Peptide_Per_Pep/Num_Cond_Peptide_Per_Pep	1	0.875	0.625
##	-logP	0.46	0.529	1.081
```

Results generated from indels. 
```
[R]
print(head(Result_HLA2_INDEL))
```

```
##	HLA	Pos	Gene	Evaluated_Mutant_Peptide_Core	Evaluated_Mutant_Peptide	Mut_IC50	Mut_Rank	Chr	NM_ID	Change	Ref	Alt	Prob	Mutation_Prob.	Exon_Start	Exon_End	Mutation_Position	Total_Depth	Tumor_Depth	Wt_Peptide	Mutant_Peptide	Total_RNA	Tumor_RNA_Ratio	Tumor_RNA	Tumor_RNA_based_on_DNA	MutRatio	MutRatio_Min	MutRatio_Max
##	HLA-DPA10103-DPB10201	3	0_UGT2B28	YHGIPMVGI	EAIYHGIPMVGIPLV	682.35	28	4	NM_053039	Out_c.1186dupT	-	T	0	0	70146192	70160768	4_70156404	84	43	EAIYHGIPMVGIPLFWDQPDNIAHMKAKGA	EAIYHGIPMVGIPLVLGSTX	0	NA	NA	0	NA	NA	NA
##	HLA-DPA10103-DPB10201	2	0_UGT2B28	YHGIPMVGI	AIYHGIPMVGIPLVL	643.18	27	4	NM_053039	Out_c.1186dupT	-	T	0	0	70146192	70160768	4_70156404	84	43	EAIYHGIPMVGIPLFWDQPDNIAHMKAKGA	EAIYHGIPMVGIPLVLGSTX	0	NA	NA	0	NA	NA	NA
##	HLA-DPA10103-DPB10201	1	0_UGT2B28	YHGIPMVGI	IYHGIPMVGIPLVLG	799.71	31	4	NM_053039	Out_c.1186dupT	-	T	0	0	70146192	70160768	4_70156404	84	43	EAIYHGIPMVGIPLFWDQPDNIAHMKAKGA	EAIYHGIPMVGIPLVLGSTX	0	NA	NA	0	NA	NA	NA
##	HLA-DPA10103-DPB10201	4	0_UGT2B28	PMVGIPLVL	YHGIPMVGIPLVLGS	984.87	35	4	NM_053039	Out_c.1186dupT	-	T	0	0	70146192	70160768	4_70156404	84	43	EAIYHGIPMVGIPLFWDQPDNIAHMKAKGA	EAIYHGIPMVGIPLVLGSTX	0	NA	NA	0	NA	NA	NA
##	HLA-DPA10103-DPB10201	3	0_UGT2B28	PMVGIPLVL	HGIPMVGIPLVLGST	1262.77	40	4	NM_053039	Out_c.1186dupT	-	T	0	0	70146192	70160768	4_70156404	84	43	EAIYHGIPMVGIPLFWDQPDNIAHMKAKGA	EAIYHGIPMVGIPLVLGSTX	0	NA	NA	0	NA	NA	NA
##	HLA-DPA10103-DPB10201	2	1_COL12A1	YQRDELLAA	QYYQRDELLAAIKKF	528.36	23	6	NM_004370	Out_c.628delA	T	-	0	0	75794041	75915769	6_75899298	195	122	QYYQRDELLAAIKKIPYKGGNTMTGDAIDYLVK	QYYQRDELLAAIKKFHIKVATQX	0	NA	NA	0	NA	NA	NA
```

Count the number of neoantigens of which threshold is IC 50 of mutatnt peptides < 500. 
```
[R]
print(Export_Summary_IndelSV(Input = Result_HLA2_INDEL, Mut_IC50_th = 500))
```

```
##	Num_All_Alteration	Num_Evaluated_Alteration	Num_Alteration_Generating_NeoAg	Num_All_Peptide	Num_Evaluated_Peptide	Num_Peptide_Generating_NeoAg
##	3	3	3	45	45	32
```

Count the number of neoantigens for each peptide of which threshold is IC 50 of mutatnt peptides < 500. 
```
[R]
print(Export_Summary_IndelSV_perFragments(Input = Result_HLA2_INDEL, Mut_IC50_th = 500))
```

```
##		EAIYHGIPMVGIPLVLGSTX-0_UGT2B28	QYYQRDELLAAIKKFHIKVATQX-1_COL12A1	IMEMARDFLPSLKNPFWKPSILPIFMYKHCSVQFSVRHGDVQTKVHX-2_SLCO1C1
##	Num_Peptide_Per_Pep	5	8	32
##	Num_Cond_Peptide_Per_Pep	5	8	32
##	Num_Rest_Peptide_Per_Pep	5	7	20
##	Num_Rest_Peptide_Per_Pep/Num_Cond_Peptide_Per_Pep	1	0.875	0.625
##	-logP	0.46	0.529	1.081
```

Results generated from SVs. 
```
[R]
print(head(Result_HLA1_SV))
```

```
##	HLA	Pos	Gene	Evaluated_Mutant_Peptide_Core	Evaluated_Mutant_Peptide	Mut_IC50	Mut_Rank	Chr	NM_ID	Change	Ref	Alt	Prob	Mutation_Prob.		Exon_Start	Exon_End	Mutation_Position	Total_Depth	Tumor_Depth	Wt_Peptide	Mutant_Peptide	Total_RNA	Tumor_RNA_Ratio	Tumor_RNA	Tumor_RNA_based_on_DNA	MutRatio	MutRatio_Min
##	HLA-A*02:01	1	0_AAR2	YNSWEVGPKF	DYNSWEVGPKFRE	41548.2	77.1323	20	NM_015511_NM_015906	In_AAR2_exon_TRIM33_exon	G	G]1:115005805]	0	0	34824338	34844863	20_34827929	0	0	MAAVQMDPELAKRLFFEGATVVILNMPKGTEFGIDYNSWEVGPKFRGVKMIPPGIHFLHYSSVDKANPKEVGPRMGFFLSLHQRGLTVLRWSTLREEVDLSPAPESEVEAMRANLQELDQFLGPYPYATLKKWISLTNFISEATVEKLQPENRQICAFSDVLPVLSMKHTKDRVGQNLPRCGIECKSYQEGLARLPEMKPRAGTEIRFSELPTQMFPEGATPAEITKHSMDLSYALETVLNKQFPSSPQDVLGELQFAFvcfLLGNVYEAFEHWKRLLNLLCRSEAAMMKHHTLYINLISILYHQLGEIPADFFVDIVSQDNFLTSTLQVFFSSACSIAVDATLRKKAEKFQAHLTKKFRWDFAAEPEDCAPVVVELPEGIEMGXXMAENKGGGEAESGGGGSGSAPVTAGAAGPAAQEAEPPLTAVLVEEEEEEGGRAGAEGGAAGPDDGGVAAASSGSAQAASSPAASVGTGVAGGAVSTPAPAPASAPAPGPSAGPPPGPPASLLDTCAVCQQSLQSRREAEPKLLPCLHSFCLRCLPEPERQLSVPIPGGSNGDIQQVGVIRCPVCRQECRQIDLVDNYFVKDTSEAPSSSDEKSEQVCTSCEDNASAVGFCVECGEWLCKTCIEAHQRVKFTKDHLIRKKEDVSESVGASGQRPVFCPVHKQEQLKLFCETCDRLTCRDCQLLEHKEHRYQFLEEAFQNQKGAIENLLAKLLEKKNYVHFAATQVQNRIKEVNETNKRVEQEIKVAIFTLINEINKKGKSLLQQLENVTKERQMKLLQQQNDITGLSRQVKHVMNFTNWAIASGSSTALLYSKRLITFQLRHILKARCDPVPAANGAIRFHCDPTFWAKNVVNLGNLVIESKPAPGYTPNVVVGQVPPGTNHISKTPGQINLAQLRLQHMQQQVYAQKHQQLQQMRMQQPPAPVPTTTTTTQQHPRQAAPQMLQQQPPRLISVQTMQRGNMNCGAFQAHQMRLAQNAARIPGIPRHSGPQYSMMQPHLQRQHSNPGHAGPFPVVSVHNTTINPTSPTTATMANANRGPTSPSVTAIELIPSVTNPENLPSLPDIPPIQLEDAGSSSLDNLLSRYISGSHLPPQPTSTMNPSPGPSALSPGSSGLSNSHTPVRPPSTSSTGSRGSCGSSGRTAEKTSLSFKSDQVKVKQEPGTEDEICSFSGGVKQEKTEDGRRSACMLSSPESSLTPPLSTNLHLESELDALASLENHVKIEPADMNESCKQSGLSSLVNGKSPIRSLMHRSARIGGDGNNKDDDPNEDWCAVCQNGGDLLCCEKCPKVFHLTCHVPTLLSFPSGDWICTFCRDIGKPEVEYDCDNLQHSKKGKTAQGLSPVDQRKCERLLLYLYCHELSIEFQEPVPASIPNYYKIIKKPMDLSTVKKKLQKKHSQHYQIPDDFVADVRLIFKNCERFNEMMKVVQVYADTQEINLKADSEVAQAGKAVALYFEDKLTEIYSDRTFAPLPEFEQEEDDGEVTEDSDEDFIQPRRKRLKSDERPVHIKX	DYNSWEVGPKFREQLKLFCETCDRLTCRDCQLLE	3.03498	NA	NA	NaN	NA	NA	NA
##	HLA-A*02:01	2	0_AAR2	YNSWEVGPKF	YNSWEVGPKFREQ	37159.3	60.6291	20	NM_015511_NM_015906	In_AAR2_exon_TRIM33_exon	G	G]1:115005805]	0	0	34824338	34844863	20_34827929	0	0	MAAVQMDPELAKRLFFEGATVVILNMPKGTEFGIDYNSWEVGPKFRGVKMIPPGIHFLHYSSVDKANPKEVGPRMGFFLSLHQRGLTVLRWSTLREEVDLSPAPESEVEAMRANLQELDQFLGPYPYATLKKWISLTNFISEATVEKLQPENRQICAFSDVLPVLSMKHTKDRVGQNLPRCGIECKSYQEGLARLPEMKPRAGTEIRFSELPTQMFPEGATPAEITKHSMDLSYALETVLNKQFPSSPQDVLGELQFAFvcfLLGNVYEAFEHWKRLLNLLCRSEAAMMKHHTLYINLISILYHQLGEIPADFFVDIVSQDNFLTSTLQVFFSSACSIAVDATLRKKAEKFQAHLTKKFRWDFAAEPEDCAPVVVELPEGIEMGXXMAENKGGGEAESGGGGSGSAPVTAGAAGPAAQEAEPPLTAVLVEEEEEEGGRAGAEGGAAGPDDGGVAAASSGSAQAASSPAASVGTGVAGGAVSTPAPAPASAPAPGPSAGPPPGPPASLLDTCAVCQQSLQSRREAEPKLLPCLHSFCLRCLPEPERQLSVPIPGGSNGDIQQVGVIRCPVCRQECRQIDLVDNYFVKDTSEAPSSSDEKSEQVCTSCEDNASAVGFCVECGEWLCKTCIEAHQRVKFTKDHLIRKKEDVSESVGASGQRPVFCPVHKQEQLKLFCETCDRLTCRDCQLLEHKEHRYQFLEEAFQNQKGAIENLLAKLLEKKNYVHFAATQVQNRIKEVNETNKRVEQEIKVAIFTLINEINKKGKSLLQQLENVTKERQMKLLQQQNDITGLSRQVKHVMNFTNWAIASGSSTALLYSKRLITFQLRHILKARCDPVPAANGAIRFHCDPTFWAKNVVNLGNLVIESKPAPGYTPNVVVGQVPPGTNHISKTPGQINLAQLRLQHMQQQVYAQKHQQLQQMRMQQPPAPVPTTTTTTQQHPRQAAPQMLQQQPPRLISVQTMQRGNMNCGAFQAHQMRLAQNAARIPGIPRHSGPQYSMMQPHLQRQHSNPGHAGPFPVVSVHNTTINPTSPTTATMANANRGPTSPSVTAIELIPSVTNPENLPSLPDIPPIQLEDAGSSSLDNLLSRYISGSHLPPQPTSTMNPSPGPSALSPGSSGLSNSHTPVRPPSTSSTGSRGSCGSSGRTAEKTSLSFKSDQVKVKQEPGTEDEICSFSGGVKQEKTEDGRRSACMLSSPESSLTPPLSTNLHLESELDALASLENHVKIEPADMNESCKQSGLSSLVNGKSPIRSLMHRSARIGGDGNNKDDDPNEDWCAVCQNGGDLLCCEKCPKVFHLTCHVPTLLSFPSGDWICTFCRDIGKPEVEYDCDNLQHSKKGKTAQGLSPVDQRKCERLLLYLYCHELSIEFQEPVPASIPNYYKIIKKPMDLSTVKKKLQKKHSQHYQIPDDFVADVRLIFKNCERFNEMMKVVQVYADTQEINLKADSEVAQAGKAVALYFEDKLTEIYSDRTFAPLPEFEQEEDDGEVTEDSDEDFIQPRRKRLKSDERPVHIKX	DYNSWEVGPKFREQLKLFCETCDRLTCRDCQLLE	3.03498	NA	NA	NaN	NA	NA	NA
##	HLA-A*02:01	3	0_AAR2	NSWEVGPKFREQL	NSWEVGPKFREQL	17541.5	22.2465	20	NM_015511_NM_015906	In_AAR2_exon_TRIM33_exon	G	G]1:115005805]	0	0	34824338	34844863	20_34827929	0	0	MAAVQMDPELAKRLFFEGATVVILNMPKGTEFGIDYNSWEVGPKFRGVKMIPPGIHFLHYSSVDKANPKEVGPRMGFFLSLHQRGLTVLRWSTLREEVDLSPAPESEVEAMRANLQELDQFLGPYPYATLKKWISLTNFISEATVEKLQPENRQICAFSDVLPVLSMKHTKDRVGQNLPRCGIECKSYQEGLARLPEMKPRAGTEIRFSELPTQMFPEGATPAEITKHSMDLSYALETVLNKQFPSSPQDVLGELQFAFvcfLLGNVYEAFEHWKRLLNLLCRSEAAMMKHHTLYINLISILYHQLGEIPADFFVDIVSQDNFLTSTLQVFFSSACSIAVDATLRKKAEKFQAHLTKKFRWDFAAEPEDCAPVVVELPEGIEMGXXMAENKGGGEAESGGGGSGSAPVTAGAAGPAAQEAEPPLTAVLVEEEEEEGGRAGAEGGAAGPDDGGVAAASSGSAQAASSPAASVGTGVAGGAVSTPAPAPASAPAPGPSAGPPPGPPASLLDTCAVCQQSLQSRREAEPKLLPCLHSFCLRCLPEPERQLSVPIPGGSNGDIQQVGVIRCPVCRQECRQIDLVDNYFVKDTSEAPSSSDEKSEQVCTSCEDNASAVGFCVECGEWLCKTCIEAHQRVKFTKDHLIRKKEDVSESVGASGQRPVFCPVHKQEQLKLFCETCDRLTCRDCQLLEHKEHRYQFLEEAFQNQKGAIENLLAKLLEKKNYVHFAATQVQNRIKEVNETNKRVEQEIKVAIFTLINEINKKGKSLLQQLENVTKERQMKLLQQQNDITGLSRQVKHVMNFTNWAIASGSSTALLYSKRLITFQLRHILKARCDPVPAANGAIRFHCDPTFWAKNVVNLGNLVIESKPAPGYTPNVVVGQVPPGTNHISKTPGQINLAQLRLQHMQQQVYAQKHQQLQQMRMQQPPAPVPTTTTTTQQHPRQAAPQMLQQQPPRLISVQTMQRGNMNCGAFQAHQMRLAQNAARIPGIPRHSGPQYSMMQPHLQRQHSNPGHAGPFPVVSVHNTTINPTSPTTATMANANRGPTSPSVTAIELIPSVTNPENLPSLPDIPPIQLEDAGSSSLDNLLSRYISGSHLPPQPTSTMNPSPGPSALSPGSSGLSNSHTPVRPPSTSSTGSRGSCGSSGRTAEKTSLSFKSDQVKVKQEPGTEDEICSFSGGVKQEKTEDGRRSACMLSSPESSLTPPLSTNLHLESELDALASLENHVKIEPADMNESCKQSGLSSLVNGKSPIRSLMHRSARIGGDGNNKDDDPNEDWCAVCQNGGDLLCCEKCPKVFHLTCHVPTLLSFPSGDWICTFCRDIGKPEVEYDCDNLQHSKKGKTAQGLSPVDQRKCERLLLYLYCHELSIEFQEPVPASIPNYYKIIKKPMDLSTVKKKLQKKHSQHYQIPDDFVADVRLIFKNCERFNEMMKVVQVYADTQEINLKADSEVAQAGKAVALYFEDKLTEIYSDRTFAPLPEFEQEEDDGEVTEDSDEDFIQPRRKRLKSDERPVHIKX	DYNSWEVGPKFREQLKLFCETCDRLTCRDCQLLE	3.03498	NA	NA	NaN	NA	NA	NA
##	HLA-A*02:01	4	0_AAR2	WEVGPKFREQLK	SWEVGPKFREQLK	41225.3	75.8035	20	NM_015511_NM_015906	In_AAR2_exon_TRIM33_exon	G	G]1:115005805]	0	0	34824338	34844863	20_34827929	0	0	MAAVQMDPELAKRLFFEGATVVILNMPKGTEFGIDYNSWEVGPKFRGVKMIPPGIHFLHYSSVDKANPKEVGPRMGFFLSLHQRGLTVLRWSTLREEVDLSPAPESEVEAMRANLQELDQFLGPYPYATLKKWISLTNFISEATVEKLQPENRQICAFSDVLPVLSMKHTKDRVGQNLPRCGIECKSYQEGLARLPEMKPRAGTEIRFSELPTQMFPEGATPAEITKHSMDLSYALETVLNKQFPSSPQDVLGELQFAFvcfLLGNVYEAFEHWKRLLNLLCRSEAAMMKHHTLYINLISILYHQLGEIPADFFVDIVSQDNFLTSTLQVFFSSACSIAVDATLRKKAEKFQAHLTKKFRWDFAAEPEDCAPVVVELPEGIEMGXXMAENKGGGEAESGGGGSGSAPVTAGAAGPAAQEAEPPLTAVLVEEEEEEGGRAGAEGGAAGPDDGGVAAASSGSAQAASSPAASVGTGVAGGAVSTPAPAPASAPAPGPSAGPPPGPPASLLDTCAVCQQSLQSRREAEPKLLPCLHSFCLRCLPEPERQLSVPIPGGSNGDIQQVGVIRCPVCRQECRQIDLVDNYFVKDTSEAPSSSDEKSEQVCTSCEDNASAVGFCVECGEWLCKTCIEAHQRVKFTKDHLIRKKEDVSESVGASGQRPVFCPVHKQEQLKLFCETCDRLTCRDCQLLEHKEHRYQFLEEAFQNQKGAIENLLAKLLEKKNYVHFAATQVQNRIKEVNETNKRVEQEIKVAIFTLINEINKKGKSLLQQLENVTKERQMKLLQQQNDITGLSRQVKHVMNFTNWAIASGSSTALLYSKRLITFQLRHILKARCDPVPAANGAIRFHCDPTFWAKNVVNLGNLVIESKPAPGYTPNVVVGQVPPGTNHISKTPGQINLAQLRLQHMQQQVYAQKHQQLQQMRMQQPPAPVPTTTTTTQQHPRQAAPQMLQQQPPRLISVQTMQRGNMNCGAFQAHQMRLAQNAARIPGIPRHSGPQYSMMQPHLQRQHSNPGHAGPFPVVSVHNTTINPTSPTTATMANANRGPTSPSVTAIELIPSVTNPENLPSLPDIPPIQLEDAGSSSLDNLLSRYISGSHLPPQPTSTMNPSPGPSALSPGSSGLSNSHTPVRPPSTSSTGSRGSCGSSGRTAEKTSLSFKSDQVKVKQEPGTEDEICSFSGGVKQEKTEDGRRSACMLSSPESSLTPPLSTNLHLESELDALASLENHVKIEPADMNESCKQSGLSSLVNGKSPIRSLMHRSARIGGDGNNKDDDPNEDWCAVCQNGGDLLCCEKCPKVFHLTCHVPTLLSFPSGDWICTFCRDIGKPEVEYDCDNLQHSKKGKTAQGLSPVDQRKCERLLLYLYCHELSIEFQEPVPASIPNYYKIIKKPMDLSTVKKKLQKKHSQHYQIPDDFVADVRLIFKNCERFNEMMKVVQVYADTQEINLKADSEVAQAGKAVALYFEDKLTEIYSDRTFAPLPEFEQEEDDGEVTEDSDEDFIQPRRKRLKSDERPVHIKX	DYNSWEVGPKFREQLKLFCETCDRLTCRDCQLLE	3.03498	NA	NA	NaN	NA	NA	NA
##	HLA-A*02:01	5	0_AAR2	WEVGPKFREQLKL	WEVGPKFREQLKL	30625.8	42.6216	20	NM_015511_NM_015906	In_AAR2_exon_TRIM33_exon	G	G]1:115005805]	0	0	34824338	34844863	20_34827929	0	0	MAAVQMDPELAKRLFFEGATVVILNMPKGTEFGIDYNSWEVGPKFRGVKMIPPGIHFLHYSSVDKANPKEVGPRMGFFLSLHQRGLTVLRWSTLREEVDLSPAPESEVEAMRANLQELDQFLGPYPYATLKKWISLTNFISEATVEKLQPENRQICAFSDVLPVLSMKHTKDRVGQNLPRCGIECKSYQEGLARLPEMKPRAGTEIRFSELPTQMFPEGATPAEITKHSMDLSYALETVLNKQFPSSPQDVLGELQFAFvcfLLGNVYEAFEHWKRLLNLLCRSEAAMMKHHTLYINLISILYHQLGEIPADFFVDIVSQDNFLTSTLQVFFSSACSIAVDATLRKKAEKFQAHLTKKFRWDFAAEPEDCAPVVVELPEGIEMGXXMAENKGGGEAESGGGGSGSAPVTAGAAGPAAQEAEPPLTAVLVEEEEEEGGRAGAEGGAAGPDDGGVAAASSGSAQAASSPAASVGTGVAGGAVSTPAPAPASAPAPGPSAGPPPGPPASLLDTCAVCQQSLQSRREAEPKLLPCLHSFCLRCLPEPERQLSVPIPGGSNGDIQQVGVIRCPVCRQECRQIDLVDNYFVKDTSEAPSSSDEKSEQVCTSCEDNASAVGFCVECGEWLCKTCIEAHQRVKFTKDHLIRKKEDVSESVGASGQRPVFCPVHKQEQLKLFCETCDRLTCRDCQLLEHKEHRYQFLEEAFQNQKGAIENLLAKLLEKKNYVHFAATQVQNRIKEVNETNKRVEQEIKVAIFTLINEINKKGKSLLQQLENVTKERQMKLLQQQNDITGLSRQVKHVMNFTNWAIASGSSTALLYSKRLITFQLRHILKARCDPVPAANGAIRFHCDPTFWAKNVVNLGNLVIESKPAPGYTPNVVVGQVPPGTNHISKTPGQINLAQLRLQHMQQQVYAQKHQQLQQMRMQQPPAPVPTTTTTTQQHPRQAAPQMLQQQPPRLISVQTMQRGNMNCGAFQAHQMRLAQNAARIPGIPRHSGPQYSMMQPHLQRQHSNPGHAGPFPVVSVHNTTINPTSPTTATMANANRGPTSPSVTAIELIPSVTNPENLPSLPDIPPIQLEDAGSSSLDNLLSRYISGSHLPPQPTSTMNPSPGPSALSPGSSGLSNSHTPVRPPSTSSTGSRGSCGSSGRTAEKTSLSFKSDQVKVKQEPGTEDEICSFSGGVKQEKTEDGRRSACMLSSPESSLTPPLSTNLHLESELDALASLENHVKIEPADMNESCKQSGLSSLVNGKSPIRSLMHRSARIGGDGNNKDDDPNEDWCAVCQNGGDLLCCEKCPKVFHLTCHVPTLLSFPSGDWICTFCRDIGKPEVEYDCDNLQHSKKGKTAQGLSPVDQRKCERLLLYLYCHELSIEFQEPVPASIPNYYKIIKKPMDLSTVKKKLQKKHSQHYQIPDDFVADVRLIFKNCERFNEMMKVVQVYADTQEINLKADSEVAQAGKAVALYFEDKLTEIYSDRTFAPLPEFEQEEDDGEVTEDSDEDFIQPRRKRLKSDERPVHIKX	DYNSWEVGPKFREQLKLFCETCDRLTCRDCQLLE	3.03498	NA	NA	NaN	NA	NA	NA
##	HLA-A*02:01	6	0_AAR2	EVGPKFREQLKL	EVGPKFREQLKLF	40083	71.2831	20	NM_015511_NM_015906	In_AAR2_exon_TRIM33_exon	G	G]1:115005805]	0	0	34824338	34844863	20_34827929	0	0	MAAVQMDPELAKRLFFEGATVVILNMPKGTEFGIDYNSWEVGPKFRGVKMIPPGIHFLHYSSVDKANPKEVGPRMGFFLSLHQRGLTVLRWSTLREEVDLSPAPESEVEAMRANLQELDQFLGPYPYATLKKWISLTNFISEATVEKLQPENRQICAFSDVLPVLSMKHTKDRVGQNLPRCGIECKSYQEGLARLPEMKPRAGTEIRFSELPTQMFPEGATPAEITKHSMDLSYALETVLNKQFPSSPQDVLGELQFAFvcfLLGNVYEAFEHWKRLLNLLCRSEAAMMKHHTLYINLISILYHQLGEIPADFFVDIVSQDNFLTSTLQVFFSSACSIAVDATLRKKAEKFQAHLTKKFRWDFAAEPEDCAPVVVELPEGIEMGXXMAENKGGGEAESGGGGSGSAPVTAGAAGPAAQEAEPPLTAVLVEEEEEEGGRAGAEGGAAGPDDGGVAAASSGSAQAASSPAASVGTGVAGGAVSTPAPAPASAPAPGPSAGPPPGPPASLLDTCAVCQQSLQSRREAEPKLLPCLHSFCLRCLPEPERQLSVPIPGGSNGDIQQVGVIRCPVCRQECRQIDLVDNYFVKDTSEAPSSSDEKSEQVCTSCEDNASAVGFCVECGEWLCKTCIEAHQRVKFTKDHLIRKKEDVSESVGASGQRPVFCPVHKQEQLKLFCETCDRLTCRDCQLLEHKEHRYQFLEEAFQNQKGAIENLLAKLLEKKNYVHFAATQVQNRIKEVNETNKRVEQEIKVAIFTLINEINKKGKSLLQQLENVTKERQMKLLQQQNDITGLSRQVKHVMNFTNWAIASGSSTALLYSKRLITFQLRHILKARCDPVPAANGAIRFHCDPTFWAKNVVNLGNLVIESKPAPGYTPNVVVGQVPPGTNHISKTPGQINLAQLRLQHMQQQVYAQKHQQLQQMRMQQPPAPVPTTTTTTQQHPRQAAPQMLQQQPPRLISVQTMQRGNMNCGAFQAHQMRLAQNAARIPGIPRHSGPQYSMMQPHLQRQHSNPGHAGPFPVVSVHNTTINPTSPTTATMANANRGPTSPSVTAIELIPSVTNPENLPSLPDIPPIQLEDAGSSSLDNLLSRYISGSHLPPQPTSTMNPSPGPSALSPGSSGLSNSHTPVRPPSTSSTGSRGSCGSSGRTAEKTSLSFKSDQVKVKQEPGTEDEICSFSGGVKQEKTEDGRRSACMLSSPESSLTPPLSTNLHLESELDALASLENHVKIEPADMNESCKQSGLSSLVNGKSPIRSLMHRSARIGGDGNNKDDDPNEDWCAVCQNGGDLLCCEKCPKVFHLTCHVPTLLSFPSGDWICTFCRDIGKPEVEYDCDNLQHSKKGKTAQGLSPVDQRKCERLLLYLYCHELSIEFQEPVPASIPNYYKIIKKPMDLSTVKKKLQKKHSQHYQIPDDFVADVRLIFKNCERFNEMMKVVQVYADTQEINLKADSEVAQAGKAVALYFEDKLTEIYSDRTFAPLPEFEQEEDDGEVTEDSDEDFIQPRRKRLKSDERPVHIKX	DYNSWEVGPKFREQLKLFCETCDRLTCRDCQLLE	3.03498	NA	NA	NaN	NA	NA	NA
```

Count the number of neoantigens of which threshold is IC 50 of mutatnt peptides < 500. 
```
[R]
print(Export_Summary_IndelSV(Result_HLA1_SV, Mut_IC50_th = 500))
```

```
##	Num_All_Alteration	Num_Evaluated_Alteration	Num_Alteration_Generating_NeoAg	Num_All_Peptide	Num_Evaluated_Peptide	Num_Peptide_Generating_NeoAg	
##	16	16	10	3183	3183	119
```

Count the number of neoantigens for each peptide of which threshold is IC 50 of mutatnt peptides < 500. 
```
[R]
print(Export_Summary_IndelSV_perFragments(Result_HLA1_SV, Mut_IC50_th = 500))
```

```
##		DYNSWEVGPKFREQLKLFCETCDRLTCRDCQLLE-0_AAR2	GDPHELRNIFLQLSAVQEAQLKRLEVTRPRVLGSREQGQVPRMARQPPPPWVHAAFLLCLLSLGG-1_SLC25A12	TMAEKRQLFIEMLYRLX-2_DTNB	TDWLSDCCFHPSERCRCTLYGHTDSVNSIEFFP-3_SPAG16	ACLPGEEGTAERSCKRNLSTTTGSSERX-4_TMCC1	GQLRGLIPSSPQKTLECKQQKIILLRREMKETX-9_PDLIM3	GDQYIVDMANTKRTLLWTRPMTVILGKPFCVMPKQQKTAHIGFLQHIPRLSPKPCLPKLNLMMRKQRMSQNGKNVKFEESHLRAVCMSGRGMGQVWVFFLCSX-10_WDR70	LCPASRALEEKKGHELLGGPSLGAPRPGSQERTGENTAACQDHRFWAGQTAGCGRERIPCRRRQSAYQVDGIGINFTQNLYPPEX-11_EGFR	LTPMPEGLSQQQDMLKNTSKGHPDRLPLQMALT-12_ARHGEF10	PGVQGKKVRKPSPTTSASSSSWTSAGTSWPAGRRTGTATSGTATTTSVWPGCGTRMWSTQWSSVPRSRSCCSRPATTPPSKPGAPHAPCASSRHLAHGLAPSSPGLPARGAEVCWVHWSHRDPLRTSPGSVAFSRAGEVEMLIAVTPX-13_NOTCH1	SRGRGNNNRKGREVTPL-14_UBAP2	IPSSQLAAKLLHMLTMRMLSKSATGRWX-15_RASGRP2	RKSPPEKKLRRYPPGQGATIGVDFMIKTVEINGE-16_GAB2	MSVIVLPQPDEVLNLVQSYVTLRVPLYVSYVFH-17_MTRF1	LWKCLRKPVLNDRNLQLHTDKGSFLKEKNKKLKKK-18_EVI2B	RTFLQSTDPGDTGGVFEGSTLCPSDPAFSRTDDPX-19_IKZF3	NVLMVHKRSHTGGVFEGSTLCPSDPAFSRTDDPX-20_IKZF3	DALTGHLRTHSGGVFEGSTLCPSDPAFSRTDDPX-21_IKZF3	FNVLMVHKRSHTGGVFEGSTLCPSDPAFSRTDDPX-22_IKZF3	EGPANEDEDIGGGVFEGSTLCPSDPAFSRTDDPX-23_IKZF3	GEGPANEDEDIGGGVFEGSTLCPSDPAFSRTDDPX-24_IKZF3	RDALTGHLRTHSGGVFEGSTLCPSDPAFSRTDDPX-25_IKZF3
##	Num_Peptide_Per_Pep	57	237	24	51	90	120	540	432	51	810	30	90	57	51	63	132	126	126	132	126	132	132
##	Num_Cond_Peptide_Per_Pep	57	237	24	51	90	120	540	432	51	810	30	90	57	51	63	132	126	126	132	126	132	132
##	Num_Rest_Peptide_Per_Pep	0	6	6	0	0	1	39	3	0	41	0	11	0	2	1	2	4	4	4	1	1	4
##	Num_Rest_Peptide_Per_Pep/Num_Cond_Peptide_Per_Pep	0	0.025	0.25	0	0	0.008	0.072	0.007	0	0.051	0	0.122	0	0.039	0.016	0.015	0.032	0.032	0.03	0.008	0.008	0.03
##	-logP	0.782	1.496	0.391	0.759	0.644	0.759	2.37	1.956	0.759	3.405	0.391	0.644	0.782	0.759	0.805	0.805	0.782	0.782	0.805	0.782	0.805	0.805
```

Results generated from SVs. 
```
[R]
print(head(Result_HLA2_SV))
```

```
##	HLA	Pos	Gene	Evaluated_Mutant_Peptide_Core	Evaluated_Mutant_Peptide	Mut_IC50	Mut_Rank	Chr	NM_ID	Change	Ref	Alt	Prob	Mutation_Prob.	Exon_Start	Exon_End	Mutation_Position	Total_Depth	Tumor_Depth	Wt_Peptide	Mutant_Peptide	Total_RNA	Tumor_RNA_Ratio	Tumor_RNA	Tumor_RNA_based_on_DNA	MutRatio	MutRatio_Min	MutRatio_Max
##	HLA-DPA10103-DPB10201	6	0_AAR2	WEVGPKFRE	GIDYNSWEVGPKFRE	3510.19	65	20	NM_015511_NM_015906	In_AAR2_exon_TRIM33_exon	G	G]1:115005805]	0	0	34824338	34844863	20_34827929	0	0	MAAVQMDPELAKRLFFEGATVVILNMPKGTEFGIDYNSWEVGPKFRGVKMIPPGIHFLHYSSVDKANPKEVGPRMGFFLSLHQRGLTVLRWSTLREEVDLSPAPESEVEAMRANLQELDQFLGPYPYATLKKWISLTNFISEATVEKLQPENRQICAFSDVLPVLSMKHTKDRVGQNLPRCGIECKSYQEGLARLPEMKPRAGTEIRFSELPTQMFPEGATPAEITKHSMDLSYALETVLNKQFPSSPQDVLGELQFAFvcfLLGNVYEAFEHWKRLLNLLCRSEAAMMKHHTLYINLISILYHQLGEIPADFFVDIVSQDNFLTSTLQVFFSSACSIAVDATLRKKAEKFQAHLTKKFRWDFAAEPEDCAPVVVELPEGIEMGXXMAENKGGGEAESGGGGSGSAPVTAGAAGPAAQEAEPPLTAVLVEEEEEEGGRAGAEGGAAGPDDGGVAAASSGSAQAASSPAASVGTGVAGGAVSTPAPAPASAPAPGPSAGPPPGPPASLLDTCAVCQQSLQSRREAEPKLLPCLHSFCLRCLPEPERQLSVPIPGGSNGDIQQVGVIRCPVCRQECRQIDLVDNYFVKDTSEAPSSSDEKSEQVCTSCEDNASAVGFCVECGEWLCKTCIEAHQRVKFTKDHLIRKKEDVSESVGASGQRPVFCPVHKQEQLKLFCETCDRLTCRDCQLLEHKEHRYQFLEEAFQNQKGAIENLLAKLLEKKNYVHFAATQVQNRIKEVNETNKRVEQEIKVAIFTLINEINKKGKSLLQQLENVTKERQMKLLQQQNDITGLSRQVKHVMNFTNWAIASGSSTALLYSKRLITFQLRHILKARCDPVPAANGAIRFHCDPTFWAKNVVNLGNLVIESKPAPGYTPNVVVGQVPPGTNHISKTPGQINLAQLRLQHMQQQVYAQKHQQLQQMRMQQPPAPVPTTTTTTQQHPRQAAPQMLQQQPPRLISVQTMQRGNMNCGAFQAHQMRLAQNAARIPGIPRHSGPQYSMMQPHLQRQHSNPGHAGPFPVVSVHNTTINPTSPTTATMANANRGPTSPSVTAIELIPSVTNPENLPSLPDIPPIQLEDAGSSSLDNLLSRYISGSHLPPQPTSTMNPSPGPSALSPGSSGLSNSHTPVRPPSTSSTGSRGSCGSSGRTAEKTSLSFKSDQVKVKQEPGTEDEICSFSGGVKQEKTEDGRRSACMLSSPESSLTPPLSTNLHLESELDALASLENHVKIEPADMNESCKQSGLSSLVNGKSPIRSLMHRSARIGGDGNNKDDDPNEDWCAVCQNGGDLLCCEKCPKVFHLTCHVPTLLSFPSGDWICTFCRDIGKPEVEYDCDNLQHSKKGKTAQGLSPVDQRKCERLLLYLYCHELSIEFQEPVPASIPNYYKIIKKPMDLSTVKKKLQKKHSQHYQIPDDFVADVRLIFKNCERFNEMMKVVQVYADTQEINLKADSEVAQAGKAVALYFEDKLTEIYSDRTFAPLPEFEQEEDDGEVTEDSDEDFIQPRRKRLKSDERPVHIKX	GIDYNSWEVGPKFREQLKLFCETCDRLTCRDCQLLEHK	3.03498	NA	NA	NaN	NA	NA	NA
##	HLA-DPA10103-DPB10201	5	0_AAR2	WEVGPKFRE	IDYNSWEVGPKFREQ	3548.49	65	20	NM_015511_NM_015906	In_AAR2_exon_TRIM33_exon	G	G]1:115005805]	0	0	34824338	34844863	20_34827929	0	0	MAAVQMDPELAKRLFFEGATVVILNMPKGTEFGIDYNSWEVGPKFRGVKMIPPGIHFLHYSSVDKANPKEVGPRMGFFLSLHQRGLTVLRWSTLREEVDLSPAPESEVEAMRANLQELDQFLGPYPYATLKKWISLTNFISEATVEKLQPENRQICAFSDVLPVLSMKHTKDRVGQNLPRCGIECKSYQEGLARLPEMKPRAGTEIRFSELPTQMFPEGATPAEITKHSMDLSYALETVLNKQFPSSPQDVLGELQFAFvcfLLGNVYEAFEHWKRLLNLLCRSEAAMMKHHTLYINLISILYHQLGEIPADFFVDIVSQDNFLTSTLQVFFSSACSIAVDATLRKKAEKFQAHLTKKFRWDFAAEPEDCAPVVVELPEGIEMGXXMAENKGGGEAESGGGGSGSAPVTAGAAGPAAQEAEPPLTAVLVEEEEEEGGRAGAEGGAAGPDDGGVAAASSGSAQAASSPAASVGTGVAGGAVSTPAPAPASAPAPGPSAGPPPGPPASLLDTCAVCQQSLQSRREAEPKLLPCLHSFCLRCLPEPERQLSVPIPGGSNGDIQQVGVIRCPVCRQECRQIDLVDNYFVKDTSEAPSSSDEKSEQVCTSCEDNASAVGFCVECGEWLCKTCIEAHQRVKFTKDHLIRKKEDVSESVGASGQRPVFCPVHKQEQLKLFCETCDRLTCRDCQLLEHKEHRYQFLEEAFQNQKGAIENLLAKLLEKKNYVHFAATQVQNRIKEVNETNKRVEQEIKVAIFTLINEINKKGKSLLQQLENVTKERQMKLLQQQNDITGLSRQVKHVMNFTNWAIASGSSTALLYSKRLITFQLRHILKARCDPVPAANGAIRFHCDPTFWAKNVVNLGNLVIESKPAPGYTPNVVVGQVPPGTNHISKTPGQINLAQLRLQHMQQQVYAQKHQQLQQMRMQQPPAPVPTTTTTTQQHPRQAAPQMLQQQPPRLISVQTMQRGNMNCGAFQAHQMRLAQNAARIPGIPRHSGPQYSMMQPHLQRQHSNPGHAGPFPVVSVHNTTINPTSPTTATMANANRGPTSPSVTAIELIPSVTNPENLPSLPDIPPIQLEDAGSSSLDNLLSRYISGSHLPPQPTSTMNPSPGPSALSPGSSGLSNSHTPVRPPSTSSTGSRGSCGSSGRTAEKTSLSFKSDQVKVKQEPGTEDEICSFSGGVKQEKTEDGRRSACMLSSPESSLTPPLSTNLHLESELDALASLENHVKIEPADMNESCKQSGLSSLVNGKSPIRSLMHRSARIGGDGNNKDDDPNEDWCAVCQNGGDLLCCEKCPKVFHLTCHVPTLLSFPSGDWICTFCRDIGKPEVEYDCDNLQHSKKGKTAQGLSPVDQRKCERLLLYLYCHELSIEFQEPVPASIPNYYKIIKKPMDLSTVKKKLQKKHSQHYQIPDDFVADVRLIFKNCERFNEMMKVVQVYADTQEINLKADSEVAQAGKAVALYFEDKLTEIYSDRTFAPLPEFEQEEDDGEVTEDSDEDFIQPRRKRLKSDERPVHIKX	GIDYNSWEVGPKFREQLKLFCETCDRLTCRDCQLLEHK	3.03498	NA	NA	NaN	NA	NA	NA
##	HLA-DPA10103-DPB10201	4	0_AAR2	WEVGPKFRE	DYNSWEVGPKFREQL	3382.78	65	20	NM_015511_NM_015906	In_AAR2_exon_TRIM33_exon	G	G]1:115005805]	0	0	34824338	34844863	20_34827929	0	0	MAAVQMDPELAKRLFFEGATVVILNMPKGTEFGIDYNSWEVGPKFRGVKMIPPGIHFLHYSSVDKANPKEVGPRMGFFLSLHQRGLTVLRWSTLREEVDLSPAPESEVEAMRANLQELDQFLGPYPYATLKKWISLTNFISEATVEKLQPENRQICAFSDVLPVLSMKHTKDRVGQNLPRCGIECKSYQEGLARLPEMKPRAGTEIRFSELPTQMFPEGATPAEITKHSMDLSYALETVLNKQFPSSPQDVLGELQFAFvcfLLGNVYEAFEHWKRLLNLLCRSEAAMMKHHTLYINLISILYHQLGEIPADFFVDIVSQDNFLTSTLQVFFSSACSIAVDATLRKKAEKFQAHLTKKFRWDFAAEPEDCAPVVVELPEGIEMGXXMAENKGGGEAESGGGGSGSAPVTAGAAGPAAQEAEPPLTAVLVEEEEEEGGRAGAEGGAAGPDDGGVAAASSGSAQAASSPAASVGTGVAGGAVSTPAPAPASAPAPGPSAGPPPGPPASLLDTCAVCQQSLQSRREAEPKLLPCLHSFCLRCLPEPERQLSVPIPGGSNGDIQQVGVIRCPVCRQECRQIDLVDNYFVKDTSEAPSSSDEKSEQVCTSCEDNASAVGFCVECGEWLCKTCIEAHQRVKFTKDHLIRKKEDVSESVGASGQRPVFCPVHKQEQLKLFCETCDRLTCRDCQLLEHKEHRYQFLEEAFQNQKGAIENLLAKLLEKKNYVHFAATQVQNRIKEVNETNKRVEQEIKVAIFTLINEINKKGKSLLQQLENVTKERQMKLLQQQNDITGLSRQVKHVMNFTNWAIASGSSTALLYSKRLITFQLRHILKARCDPVPAANGAIRFHCDPTFWAKNVVNLGNLVIESKPAPGYTPNVVVGQVPPGTNHISKTPGQINLAQLRLQHMQQQVYAQKHQQLQQMRMQQPPAPVPTTTTTTQQHPRQAAPQMLQQQPPRLISVQTMQRGNMNCGAFQAHQMRLAQNAARIPGIPRHSGPQYSMMQPHLQRQHSNPGHAGPFPVVSVHNTTINPTSPTTATMANANRGPTSPSVTAIELIPSVTNPENLPSLPDIPPIQLEDAGSSSLDNLLSRYISGSHLPPQPTSTMNPSPGPSALSPGSSGLSNSHTPVRPPSTSSTGSRGSCGSSGRTAEKTSLSFKSDQVKVKQEPGTEDEICSFSGGVKQEKTEDGRRSACMLSSPESSLTPPLSTNLHLESELDALASLENHVKIEPADMNESCKQSGLSSLVNGKSPIRSLMHRSARIGGDGNNKDDDPNEDWCAVCQNGGDLLCCEKCPKVFHLTCHVPTLLSFPSGDWICTFCRDIGKPEVEYDCDNLQHSKKGKTAQGLSPVDQRKCERLLLYLYCHELSIEFQEPVPASIPNYYKIIKKPMDLSTVKKKLQKKHSQHYQIPDDFVADVRLIFKNCERFNEMMKVVQVYADTQEINLKADSEVAQAGKAVALYFEDKLTEIYSDRTFAPLPEFEQEEDDGEVTEDSDEDFIQPRRKRLKSDERPVHIKX	GIDYNSWEVGPKFREQLKLFCETCDRLTCRDCQLLEHK	3.03498	NA	NA	NaN	NA	NA	NA
##	HLA-DPA10103-DPB10201	3	0_AAR2	WEVGPKFRE	YNSWEVGPKFREQLK	3259.68	65	20	NM_015511_NM_015906	In_AAR2_exon_TRIM33_exon	G	G]1:115005805]	0	0	34824338	34844863	20_34827929	0	0	MAAVQMDPELAKRLFFEGATVVILNMPKGTEFGIDYNSWEVGPKFRGVKMIPPGIHFLHYSSVDKANPKEVGPRMGFFLSLHQRGLTVLRWSTLREEVDLSPAPESEVEAMRANLQELDQFLGPYPYATLKKWISLTNFISEATVEKLQPENRQICAFSDVLPVLSMKHTKDRVGQNLPRCGIECKSYQEGLARLPEMKPRAGTEIRFSELPTQMFPEGATPAEITKHSMDLSYALETVLNKQFPSSPQDVLGELQFAFvcfLLGNVYEAFEHWKRLLNLLCRSEAAMMKHHTLYINLISILYHQLGEIPADFFVDIVSQDNFLTSTLQVFFSSACSIAVDATLRKKAEKFQAHLTKKFRWDFAAEPEDCAPVVVELPEGIEMGXXMAENKGGGEAESGGGGSGSAPVTAGAAGPAAQEAEPPLTAVLVEEEEEEGGRAGAEGGAAGPDDGGVAAASSGSAQAASSPAASVGTGVAGGAVSTPAPAPASAPAPGPSAGPPPGPPASLLDTCAVCQQSLQSRREAEPKLLPCLHSFCLRCLPEPERQLSVPIPGGSNGDIQQVGVIRCPVCRQECRQIDLVDNYFVKDTSEAPSSSDEKSEQVCTSCEDNASAVGFCVECGEWLCKTCIEAHQRVKFTKDHLIRKKEDVSESVGASGQRPVFCPVHKQEQLKLFCETCDRLTCRDCQLLEHKEHRYQFLEEAFQNQKGAIENLLAKLLEKKNYVHFAATQVQNRIKEVNETNKRVEQEIKVAIFTLINEINKKGKSLLQQLENVTKERQMKLLQQQNDITGLSRQVKHVMNFTNWAIASGSSTALLYSKRLITFQLRHILKARCDPVPAANGAIRFHCDPTFWAKNVVNLGNLVIESKPAPGYTPNVVVGQVPPGTNHISKTPGQINLAQLRLQHMQQQVYAQKHQQLQQMRMQQPPAPVPTTTTTTQQHPRQAAPQMLQQQPPRLISVQTMQRGNMNCGAFQAHQMRLAQNAARIPGIPRHSGPQYSMMQPHLQRQHSNPGHAGPFPVVSVHNTTINPTSPTTATMANANRGPTSPSVTAIELIPSVTNPENLPSLPDIPPIQLEDAGSSSLDNLLSRYISGSHLPPQPTSTMNPSPGPSALSPGSSGLSNSHTPVRPPSTSSTGSRGSCGSSGRTAEKTSLSFKSDQVKVKQEPGTEDEICSFSGGVKQEKTEDGRRSACMLSSPESSLTPPLSTNLHLESELDALASLENHVKIEPADMNESCKQSGLSSLVNGKSPIRSLMHRSARIGGDGNNKDDDPNEDWCAVCQNGGDLLCCEKCPKVFHLTCHVPTLLSFPSGDWICTFCRDIGKPEVEYDCDNLQHSKKGKTAQGLSPVDQRKCERLLLYLYCHELSIEFQEPVPASIPNYYKIIKKPMDLSTVKKKLQKKHSQHYQIPDDFVADVRLIFKNCERFNEMMKVVQVYADTQEINLKADSEVAQAGKAVALYFEDKLTEIYSDRTFAPLPEFEQEEDDGEVTEDSDEDFIQPRRKRLKSDERPVHIKX	GIDYNSWEVGPKFREQLKLFCETCDRLTCRDCQLLEHK	3.03498	NA	NA	NaN	NA	NA	NA
##	HLA-DPA10103-DPB10201	2	0_AAR2	WEVGPKFRE	NSWEVGPKFREQLKL	3578.45	65	20	NM_015511_NM_015906	In_AAR2_exon_TRIM33_exon	G	G]1:115005805]	0	0	34824338	34844863	20_34827929	0	0	MAAVQMDPELAKRLFFEGATVVILNMPKGTEFGIDYNSWEVGPKFRGVKMIPPGIHFLHYSSVDKANPKEVGPRMGFFLSLHQRGLTVLRWSTLREEVDLSPAPESEVEAMRANLQELDQFLGPYPYATLKKWISLTNFISEATVEKLQPENRQICAFSDVLPVLSMKHTKDRVGQNLPRCGIECKSYQEGLARLPEMKPRAGTEIRFSELPTQMFPEGATPAEITKHSMDLSYALETVLNKQFPSSPQDVLGELQFAFvcfLLGNVYEAFEHWKRLLNLLCRSEAAMMKHHTLYINLISILYHQLGEIPADFFVDIVSQDNFLTSTLQVFFSSACSIAVDATLRKKAEKFQAHLTKKFRWDFAAEPEDCAPVVVELPEGIEMGXXMAENKGGGEAESGGGGSGSAPVTAGAAGPAAQEAEPPLTAVLVEEEEEEGGRAGAEGGAAGPDDGGVAAASSGSAQAASSPAASVGTGVAGGAVSTPAPAPASAPAPGPSAGPPPGPPASLLDTCAVCQQSLQSRREAEPKLLPCLHSFCLRCLPEPERQLSVPIPGGSNGDIQQVGVIRCPVCRQECRQIDLVDNYFVKDTSEAPSSSDEKSEQVCTSCEDNASAVGFCVECGEWLCKTCIEAHQRVKFTKDHLIRKKEDVSESVGASGQRPVFCPVHKQEQLKLFCETCDRLTCRDCQLLEHKEHRYQFLEEAFQNQKGAIENLLAKLLEKKNYVHFAATQVQNRIKEVNETNKRVEQEIKVAIFTLINEINKKGKSLLQQLENVTKERQMKLLQQQNDITGLSRQVKHVMNFTNWAIASGSSTALLYSKRLITFQLRHILKARCDPVPAANGAIRFHCDPTFWAKNVVNLGNLVIESKPAPGYTPNVVVGQVPPGTNHISKTPGQINLAQLRLQHMQQQVYAQKHQQLQQMRMQQPPAPVPTTTTTTQQHPRQAAPQMLQQQPPRLISVQTMQRGNMNCGAFQAHQMRLAQNAARIPGIPRHSGPQYSMMQPHLQRQHSNPGHAGPFPVVSVHNTTINPTSPTTATMANANRGPTSPSVTAIELIPSVTNPENLPSLPDIPPIQLEDAGSSSLDNLLSRYISGSHLPPQPTSTMNPSPGPSALSPGSSGLSNSHTPVRPPSTSSTGSRGSCGSSGRTAEKTSLSFKSDQVKVKQEPGTEDEICSFSGGVKQEKTEDGRRSACMLSSPESSLTPPLSTNLHLESELDALASLENHVKIEPADMNESCKQSGLSSLVNGKSPIRSLMHRSARIGGDGNNKDDDPNEDWCAVCQNGGDLLCCEKCPKVFHLTCHVPTLLSFPSGDWICTFCRDIGKPEVEYDCDNLQHSKKGKTAQGLSPVDQRKCERLLLYLYCHELSIEFQEPVPASIPNYYKIIKKPMDLSTVKKKLQKKHSQHYQIPDDFVADVRLIFKNCERFNEMMKVVQVYADTQEINLKADSEVAQAGKAVALYFEDKLTEIYSDRTFAPLPEFEQEEDDGEVTEDSDEDFIQPRRKRLKSDERPVHIKX	GIDYNSWEVGPKFREQLKLFCETCDRLTCRDCQLLEHK	3.03498	NA	NA	NaN	NA	NA	NA
##	HLA-DPA10103-DPB10201	6	0_AAR2	KFREQLKLF	SWEVGPKFREQLKLF	2000.52	50	20	NM_015511_NM_015906	In_AAR2_exon_TRIM33_exon	G	G]1:115005805]	0	0	34824338	34844863	20_34827929	0	0	MAAVQMDPELAKRLFFEGATVVILNMPKGTEFGIDYNSWEVGPKFRGVKMIPPGIHFLHYSSVDKANPKEVGPRMGFFLSLHQRGLTVLRWSTLREEVDLSPAPESEVEAMRANLQELDQFLGPYPYATLKKWISLTNFISEATVEKLQPENRQICAFSDVLPVLSMKHTKDRVGQNLPRCGIECKSYQEGLARLPEMKPRAGTEIRFSELPTQMFPEGATPAEITKHSMDLSYALETVLNKQFPSSPQDVLGELQFAFvcfLLGNVYEAFEHWKRLLNLLCRSEAAMMKHHTLYINLISILYHQLGEIPADFFVDIVSQDNFLTSTLQVFFSSACSIAVDATLRKKAEKFQAHLTKKFRWDFAAEPEDCAPVVVELPEGIEMGXXMAENKGGGEAESGGGGSGSAPVTAGAAGPAAQEAEPPLTAVLVEEEEEEGGRAGAEGGAAGPDDGGVAAASSGSAQAASSPAASVGTGVAGGAVSTPAPAPASAPAPGPSAGPPPGPPASLLDTCAVCQQSLQSRREAEPKLLPCLHSFCLRCLPEPERQLSVPIPGGSNGDIQQVGVIRCPVCRQECRQIDLVDNYFVKDTSEAPSSSDEKSEQVCTSCEDNASAVGFCVECGEWLCKTCIEAHQRVKFTKDHLIRKKEDVSESVGASGQRPVFCPVHKQEQLKLFCETCDRLTCRDCQLLEHKEHRYQFLEEAFQNQKGAIENLLAKLLEKKNYVHFAATQVQNRIKEVNETNKRVEQEIKVAIFTLINEINKKGKSLLQQLENVTKERQMKLLQQQNDITGLSRQVKHVMNFTNWAIASGSSTALLYSKRLITFQLRHILKARCDPVPAANGAIRFHCDPTFWAKNVVNLGNLVIESKPAPGYTPNVVVGQVPPGTNHISKTPGQINLAQLRLQHMQQQVYAQKHQQLQQMRMQQPPAPVPTTTTTTQQHPRQAAPQMLQQQPPRLISVQTMQRGNMNCGAFQAHQMRLAQNAARIPGIPRHSGPQYSMMQPHLQRQHSNPGHAGPFPVVSVHNTTINPTSPTTATMANANRGPTSPSVTAIELIPSVTNPENLPSLPDIPPIQLEDAGSSSLDNLLSRYISGSHLPPQPTSTMNPSPGPSALSPGSSGLSNSHTPVRPPSTSSTGSRGSCGSSGRTAEKTSLSFKSDQVKVKQEPGTEDEICSFSGGVKQEKTEDGRRSACMLSSPESSLTPPLSTNLHLESELDALASLENHVKIEPADMNESCKQSGLSSLVNGKSPIRSLMHRSARIGGDGNNKDDDPNEDWCAVCQNGGDLLCCEKCPKVFHLTCHVPTLLSFPSGDWICTFCRDIGKPEVEYDCDNLQHSKKGKTAQGLSPVDQRKCERLLLYLYCHELSIEFQEPVPASIPNYYKIIKKPMDLSTVKKKLQKKHSQHYQIPDDFVADVRLIFKNCERFNEMMKVVQVYADTQEINLKADSEVAQAGKAVALYFEDKLTEIYSDRTFAPLPEFEQEEDDGEVTEDSDEDFIQPRRKRLKSDERPVHIKX	GIDYNSWEVGPKFREQLKLFCETCDRLTCRDCQLLEHK	3.03498	NA	NA	NaN	NA	NA	NA
```

Count the number of neoantigens of which threshold is IC 50 of mutatnt peptides < 500. 
```
[R]
print(Export_Summary_IndelSV(Result_HLA2_SV, Mut_IC50_th = 500))
```

```
##	Num_All_Alteration	Num_Evaluated_Alteration	Num_Alteration_Generating_NeoAg	Num_All_Peptide	Num_Evaluated_Peptide	Num_Peptide_Generating_NeoAg
##	16	16	12	588	588	298
```

Count the number of neoantigens for each peptide of which threshold is IC 50 of mutatnt peptides < 500. 
```
[R]
print(Export_Summary_IndelSV_perFragments(Result_HLA2_SV, Mut_IC50_th = 500))
```

```
##		GIDYNSWEVGPKFREQLKLFCETCDRLTCRDCQLLEHK-0_AAR2	KRGDPHELRNIFLQLSAVQEAQLKRLEVTRPRVLGSREQGQVPRMARQPPPPWVHAAFLLCLLSLGGAI-1_SLC25A12	RKTMAEKRQLFIEMLYRLX-2_DTNB	GHTDWLSDCCFHPSERCRCTLYGHTDSVNSIEFFPFS-3_SPAG16	AAACLPGEEGTAERSCKRNLSTTTGSSERX-4_TMCC1	MKSCKRNLSTTTGSSERX-5_TMCC1	MVQRFSLRRQLSKSCKRNLSTTTGSSERX-6_TMCC1	PMQRKGSYCWDFEESCKRNLSTTTGSSERX-7_TMCC1	HNIRPKPFVIPGRSRTLECKQQKIILLRREMKETX-8_PDLIM3	IKGDQYIVDMANTKRTLLWTRPMTVILGKPFCVMPKQQKTAHIGFLQHIPRLSPKPCLPKLNLMMRKQRMSQNGKNVKFEESHLRAVCMSGRGMGQVWVFFLCSX-9_WDR70	AALCPASRALEEKKGHELLGGPSLGAPRPGSQERTGENTAACQDHRFWAGQTAGCGRERIPCRRRQSAYQVDGIGINFTQNLYPPEX-10_EGFR	AILTPMPEGLSQQQDMLKNTSKGHPDRLPLQMALTEL-11_ARHGEF10	LKPGVQGKKVRKPSPTTSASSSSWTSAGTSWPAGRRTGTATSGTATTTSVWPGCGTRMWSTQWSSVPRSRSCCSRPATTPPSKPGAPHAPCASSRHLAHGLAPSSPGLPARGAEVCWVHWSHRDPLRTSPGSVAFSRAGEVEMLIAVTPX-12_NOTCH1	ESSRGRGNNNRKGREVTPL-13_UBAP2	WYIPSSQLAAKLLHMLTMRMLSKSATGRWX-14_RASGRP2	WLRKSPPEKKLRRYPPGQGATIGVDFMIKTVEINGEKV-15_GAB2	GTMSVIVLPQPDEVLNLVQSYVTLRVPLYVSYVFHSP-16_MTRF1	IVLWKCLRKPVLNDRNLQLHTDKGSFLKEKNKKLKKKNK-17_EVI2B	RCRTFLQSTDPGDTGGVFEGSTLCPSDPAFSRTDDPX-18_IKZF3	SFNVLMVHKRSHTGGVFEGSTLCPSDPAFSRTDDPX-19_IKZF3	SGEGPANEDEDIGGGVFEGSTLCPSDPAFSRTDDPX-20_IKZF3	DSGEGPANEDEDIGGGVFEGSTLCPSDPAFSRTDDPX-21_IKZF3	QRRDALTGHLRTHSGGVFEGSTLCPSDPAFSRTDDPX-22_IKZF3	ISFNVLMVHKRSHTGGVFEGSTLCPSDPAFSRTDDPX-23_IKZF3
##	Num_Peptide_Per_Pep	14	44	4	13	15	3	14	15	20	90	72	13	135	5	15	14	13	15	22	21	21	22	22	22
##	Num_Cond_Peptide_Per_Pep	14	44	4	13	15	3	14	15	20	90	72	13	135	5	15	14	13	15	22	21	21	22	22	22
##	Num_Rest_Peptide_Per_Pep	0	19	4	0	0	0	9	0	10	55	24	0	94	0	15	5	6	11	10	17	12	12	13	18
##	Num_Rest_Peptide_Per_Pep/Num_Cond_Peptide_Per_Pep	0	0.432	1	0	0	0	0.643	0	0.5	0.611	0.333	0	0.696	0	1	0.357	0.462	0.733	0.455	0.81	0.571	0.545	0.591	0.818
##	-logP	0.874	1.588	0.437	0.851	0.69	0.414	0.667	0.69	0.805	2.416	2.002	0.851	3.451	0.437	0.69	0.874	0.851	0.897	0.851	0.828	0.828	0.851	0.851	0.851```
```

Results generated from RNA sequences. 
```
[R]
print(head(Result_HLA1_Seq))
```

```
##	HLA	Pos	Gene	Evaluated_Mutant_Peptide_Core	Evaluated_Mutant_Peptide	Mut_IC50	Mut_Rank	Chr	NM_ID	ReadingFrame	SequenceNumber	Chrs	NM_IDs	GeneIDs	Exon_Starts	Exon_Ends	GroupID	NumOfPeptides	NumOfStops
##	HLA-A*02:01	1	0_atggcagaag	MAEDDPYLGRPEK	MAEDDPYLGRPEK	34529.5	52.6059	0		1	1	chr4;chr4	NM_001165412;NM_003998	NFKB1;NFKB1	103422515;103422515	103538459;103538459	0_1	27	0
##	HLA-A*02:01	2	0_atggcagaag	AEDDPYLGRPEKM	AEDDPYLGRPEKM	37395.3	61.456	0		1	1	chr4;chr4	NM_001165412;NM_003998	NFKB1;NFKB1	103422515;103422515	103538459;103538459	0_1	27	0
##	HLA-A*02:01	3	0_atggcagaag	YLGRPEKMF	EDDPYLGRPEKMF	41328	76.2274	0		1	1	chr4;chr4	NM_001165412;NM_003998	NFKB1;NFKB1	103422515;103422515	103538459;103538459	0_1	27	0
##	HLA-A*02:01	4	0_atggcagaag	YLGRPEKMFH	DDPYLGRPEKMFH	42864.6	82.5363	0		1	1	chr4;chr4	NM_001165412;NM_003998	NFKB1;NFKB1	103422515;103422515	103538459;103538459	0_1	27	0
##	HLA-A*02:01	5	0_atggcagaag	YLGRPEKMFHL	DPYLGRPEKMFHL	3345.9	7.5921	0		1	1	chr4;chr4	NM_001165412;NM_003998	NFKB1;NFKB1	103422515;103422515	103538459;103538459	0_1	27	0
##	HLA-A*02:01	6	0_atggcagaag	YLGRPEKMFHL	PYLGRPEKMFHLD	19715.7	24.819	0		1	1	chr4;chr4	NM_001165412;NM_003998	NFKB1;NFKB1	103422515;103422515	103538459;103538459	0_1	27	0
```

Count the number of neoantigens for each peptide of which threshold is IC 50 of mutatnt peptides < 500. 
```
[R]
print(Export_Summary_Fragments(Result_HLA1_Seq, Mut_IC50_th = 500))
```

```
## [[1]]
##                                                     atggcagaag
## Num_Peptide_Per_Grp                                     63.000
## Num_Cond_Peptide_Per_Grp                                63.000
## Num_Rest_Peptide_Per_Grp                                 5.000
## Num_Rest_Peptide_Per_Grp / Num_Cond_Peptide_Per_Grp      0.079
## -logP                                                    0.621
```

```
## [[2]]
##                                                         
## Num_Peptide_Per_NM                                63.000
## Num_Cond_Peptide_Per_NM                           63.000
## Num_Rest_Peptide_Per_NM                            5.000
## Num_Rest_Peptide_Per_NM / Num_Cond_Peptide_Per_NM  0.079
## -logP                                              0.621
```

```
## [[3]]
##                                                     MAEDDPYLGRPEKMFHLDPSLTHTIFN-0_atggcagaag-0_1
## Num_Peptide_Per_Pep                                                                       63.000
## Num_Cond_Peptide_Per_Pep                                                                  63.000
## Num_Rest_Peptide_Per_Pep                                                                   5.000
## Num_Rest_Peptide_Per_Pep / Num_Cond_Peptide_Per_Pep                                        0.079
## -logP                                                                                      0.621
```

Results generated from RNA sequences. 
```
[R]
print(head(Result_HLA2_Seq))
```

```
##	HLA	Pos	Gene	Evaluated_Mutant_Peptide_Core	Evaluated_Mutant_Peptide	Mut_IC50	Mut_Rank	Chr	NM_ID	ReadingFrame	SequenceNumber	Chrs	NM_IDs	GeneIDs	Exon_Starts	Exon_Ends	GroupID	NumOfPeptides	NumOfStops	Wt_Peptide	Mutant_Peptide	Total_RNA	Tumor_RNA_Ratio	Tumor_RNA	Tumor_RNA_based_on_DNA	MutRatio	MutRatio_Min	MutRatio_Max
##	HLA-DPA10103-DPB10201	5	0_atggcagaag	GRPEQMFHL	DDPYLGRPEQMFHLI	2116.65	55	0		1	1	chr4;chr4;chr19;chr4	NM_001165412;NM_003998;NM_005178;NM_001319226	NFKB1;NFKB1;BCL3;NFKB1	103422515;103422515;45251964;103423090	103538459;103538459;45263301;103538459	0_1	27	1	MAEDDPYLGRPEQMFHLDPSLTHTIFNPEVFQPQMALPTDGPYLQILEQPKQRGFRFRYVCEGPSHGGLPGASSEKNKKSYPQVKICNYVGPAKVIVQLVTNGKNIHLHAHSLVGKHCEDGICTVTAGPKDMVVGFANLGILHVTKKKVFETLEARMTEACIRGYNPGLLVHPDLAYLQAEGGGDRQLGDREKELIRQAALQQTKEMDLSVVRLMFTAFLPDSTGSFTRRLEPVVSDAIYDSKAPNASNLKIVRMDRTAGCVTGGEEIYLLCDKVQKDDIQIRFYEEEENGGVWEGFGDFSPTDVHRQFAIVFKTPKYKDINITKPASVFVQLRRKSDLETSEPKPFLYYPEIKDKEEVQRKRQKLMPNFSDSFGGGSGAGAGGGGMFGSGGGGGGTGSTGPGYSFPHYGFPTYGGITFHPGTTKSNAGMKHGTMDTESKKDPEGCDKSDDKNTVNLFGKVIETTEQDQEPSEATVGNGEVTLTYATGTKEESAGVQDNLFLEKAMQLAKRHANALFDYAVTGDVKMLLAVQRHLTAVQDENGDSVLHLAIIHLHSQLVRDLLEVTSGLISDDIINMRNDLYQTPLHLAVITKQEDVVEDLLRAGADLSLLDRLGNSVLHLAAKEGHDKVLSILLKHKKAALLLDHPNGDGLNAIHLAMMSNSLPCLLLLVAAGADVNAQEQKSGRTALHLAVEHDNISLAGCLLLEGDAHVDSTTYDGTTPLHIAAGRGSTRLAALLKAAGADPLVENFEPLYDLDDSWENAGEDEGVVPGTTPLDMATSWQVFDILNGKPYEPEFTSDDLLAQGDMKQLAEDVKLQLYKLLEIPDPDKNWATLAQKLGLGILNNAFRLSPAPSKTLMDNYEVSGGTVRELVEALRQMGYTEAIEVIQAASSPVKTTSQAHSLPLSPASTRQQIDELRDSDSVCDSGVETSFRKLSFTESLTSGASLLTLNKMPHDYGQEGPLEGKIXXMAEDDPYLGRPEQMFHLDPSLTHTIFNPEVFQPQMALPTADGPYLQILEQPKQRGFRFRYVCEGPSHGGLPGASSEKNKKSYPQVKICNYVGPAKVIVQLVTNGKNIHLHAHSLVGKHCEDGICTVTAGPKDMVVGFANLGILHVTKKKVFETLEARMTEACIRGYNPGLLVHPDLAYLQAEGGGDRQLGDREKELIRQAALQQTKEMDLSVVRLMFTAFLPDSTGSFTRRLEPVVSDAIYDSKAPNASNLKIVRMDRTAGCVTGGEEIYLLCDKVQKDDIQIRFYEEEENGGVWEGFGDFSPTDVHRQFAIVFKTPKYKDINITKPASVFVQLRRKSDLETSEPKPFLYYPEIKDKEEVQRKRQKLMPNFSDSFGGGSGAGAGGGGMFGSGGGGGGTGSTGPGYSFPHYGFPTYGGITFHPGTTKSNAGMKHGTMDTESKKDPEGCDKSDDKNTVNLFGKVIETTEQDQEPSEATVGNGEVTLTYATGTKEESAGVQDNLFLEKAMQLAKRHANALFDYAVTGDVKMLLAVQRHLTAVQDENGDSVLHLAIIHLHSQLVRDLLEVTSGLISDDIINMRNDLYQTPLHLAVITKQEDVVEDLLRAGADLSLLDRLGNSVLHLAAKEGHDKVLSILLKHKKAALLLDHPNGDGLNAIHLAMMSNSLPCLLLLVAAGADVNAQEQKSGRTALHLAVEHDNISLAGCLLLEGDAHVDSTTYDGTTPLHIAAGRGSTRLAALLKAAGADPLVENFEPLYDLDDSWENAGEDEGVVPGTTPLDMATSWQVFDILNGKPYEPEFTSDDLLAQGDMKQLAEDVKLQLYKLLEIPDPDKNWATLAQKLGLGILNNAFRLSPAPSKTLMDNYEVSGGTVRELVEALRQMGYTEAIEVIQAASSPVKTTSQAHSLPLSPASTRQQIDELRDSDSVCDSGVETSFRKLSFTESLTSGASLLTLNKMPHDYGQEGPLEGKIXXMPRCPAGAMDEGPVDLRTRPKAAGLPGAALPLRKRPLRAPSPEPAAPRGAAGLVVPLDPLRGGCDLPAVPGPPHGLARPEALYYPGALLPLYPTRAMGSPFPLVNLPTPLYPMMCPMEHPLSADIAMATRADEDGDTPLHIAVVQGNLPAVHRLVNLFQQGGRELDIYNNLRQTPLHLAVITTLPSVVRLLVTAGASPMALDRHGQTAAHLACEHRSPTCLRALLDSAAPGTLDLEARNYDGLTALHVAVNTECQETVQLLLERGADIDAVDIKSGRSPLIHAVENNSLSMVQLLLQHGANVNAQMYSGSSALHSASGRGLLPLVRTLVRSGADSSLKNCHNDTPLMVARSRRVIDILRGKATRPASTSQPDPSPDRSANTSPESSSRLSSNGLLSASPSSSPSQSPPRDPPGFPMAPPNFFLPSPSPPAFLPFAGVLRGPGRPVPPSPAPGGSXXMAEDDPYLGRPEQMFHLDPSLTHTIFNPEVFQPQMALPTDGPYLQILEQPKQRGFRFRYVCEGPSHGGLPGASSEKNKKSYPQVKICNYVGPAKVIVQLVTNGKNIHLHAHSLVGKHCEDGICTVTAGPKDMVVGFANLGILHVTKKKVFETLEARMTEACIRGYNPGLLVHPDLAYLQAEGGGDRQLGDREKELIRQAALQQTKEMDLSVVRLMFTAFLPDSTGSFTRRLEPVVSDAIYDSKAPNASNLKIVRMDRTAGCVTGGEEIYLLCDKVQKDDIQIRFYEEEENGGVWEGFGDFSPTDVHRQFAIVFKTPKYKDINITKPASVFVQLRRKSDLETSEPKPFLYYPEIKDKEEVQRKRQKLMPNFSDSFGGGSGAGAGGGGMFGSGGGGGGTGSTGPGYSFPHYGFPTYGGITFHPGTTKSNAGMKHGTMDTESKKDPEGCDKSDDKNTVNLFGKVIETTEQDQEPSEATVGNGEVTLTYATGTKEESAGVQDNLFLEKAMQLAKRHANALFDYAVTGDVKMLLAVQRHLTAVQDENGDSVLHLAIIHLHSQLVRDLLEVTSGLISDDIINMRNDLYQTPLHLAVITKQEDVVEDLLRAGADLSLLDRLGNSVLHLAAKEGHDKVLSILLKHKKAALLLDHPNGDGLNAIHLAMMSNSLPCLLLLVAAGADVNAQEQKSGRTALHLAVEHDNISLAGCLLLEGDAHVDSTTYDGTTPLHIAAGRGSTRLAALLKAAGADPLVENFEPLYDLDDSWENAGEDEGVVPGTTPLDMATSWQVFDILNGKPYEPEFTSDDLLAQGDMKQLAEDVKLQLYKLLEIPDPDKNWATLAQKLGLGILNNAFRLSPAPSKTLMDNYEVSGGTVRELVEALRQMGYTEAIEVIQAASSPVKTTSQAHSLPLSPASTRQQIDELRDSDSVCDSGVETSFRKLSFTESLTSGASLLTLNKMPHDYGQEGPLEGKIX	MAEDDPYLGRPEQMFHLILL	NA	NA	NA	NA	NA	NA	NA
##	HLA-DPA10103-DPB10201	4	0_atggcagaag	GRPEQMFHL	DPYLGRPEQMFHLIL	938.03	34	0		1	1	chr4;chr4;chr19;chr4	NM_001165412;NM_003998;NM_005178;NM_001319226	NFKB1;NFKB1;BCL3;NFKB1	103422515;103422515;45251964;103423090	103538459;103538459;45263301;103538459	0_1	27	1	MAEDDPYLGRPEQMFHLDPSLTHTIFNPEVFQPQMALPTDGPYLQILEQPKQRGFRFRYVCEGPSHGGLPGASSEKNKKSYPQVKICNYVGPAKVIVQLVTNGKNIHLHAHSLVGKHCEDGICTVTAGPKDMVVGFANLGILHVTKKKVFETLEARMTEACIRGYNPGLLVHPDLAYLQAEGGGDRQLGDREKELIRQAALQQTKEMDLSVVRLMFTAFLPDSTGSFTRRLEPVVSDAIYDSKAPNASNLKIVRMDRTAGCVTGGEEIYLLCDKVQKDDIQIRFYEEEENGGVWEGFGDFSPTDVHRQFAIVFKTPKYKDINITKPASVFVQLRRKSDLETSEPKPFLYYPEIKDKEEVQRKRQKLMPNFSDSFGGGSGAGAGGGGMFGSGGGGGGTGSTGPGYSFPHYGFPTYGGITFHPGTTKSNAGMKHGTMDTESKKDPEGCDKSDDKNTVNLFGKVIETTEQDQEPSEATVGNGEVTLTYATGTKEESAGVQDNLFLEKAMQLAKRHANALFDYAVTGDVKMLLAVQRHLTAVQDENGDSVLHLAIIHLHSQLVRDLLEVTSGLISDDIINMRNDLYQTPLHLAVITKQEDVVEDLLRAGADLSLLDRLGNSVLHLAAKEGHDKVLSILLKHKKAALLLDHPNGDGLNAIHLAMMSNSLPCLLLLVAAGADVNAQEQKSGRTALHLAVEHDNISLAGCLLLEGDAHVDSTTYDGTTPLHIAAGRGSTRLAALLKAAGADPLVENFEPLYDLDDSWENAGEDEGVVPGTTPLDMATSWQVFDILNGKPYEPEFTSDDLLAQGDMKQLAEDVKLQLYKLLEIPDPDKNWATLAQKLGLGILNNAFRLSPAPSKTLMDNYEVSGGTVRELVEALRQMGYTEAIEVIQAASSPVKTTSQAHSLPLSPASTRQQIDELRDSDSVCDSGVETSFRKLSFTESLTSGASLLTLNKMPHDYGQEGPLEGKIXXMAEDDPYLGRPEQMFHLDPSLTHTIFNPEVFQPQMALPTADGPYLQILEQPKQRGFRFRYVCEGPSHGGLPGASSEKNKKSYPQVKICNYVGPAKVIVQLVTNGKNIHLHAHSLVGKHCEDGICTVTAGPKDMVVGFANLGILHVTKKKVFETLEARMTEACIRGYNPGLLVHPDLAYLQAEGGGDRQLGDREKELIRQAALQQTKEMDLSVVRLMFTAFLPDSTGSFTRRLEPVVSDAIYDSKAPNASNLKIVRMDRTAGCVTGGEEIYLLCDKVQKDDIQIRFYEEEENGGVWEGFGDFSPTDVHRQFAIVFKTPKYKDINITKPASVFVQLRRKSDLETSEPKPFLYYPEIKDKEEVQRKRQKLMPNFSDSFGGGSGAGAGGGGMFGSGGGGGGTGSTGPGYSFPHYGFPTYGGITFHPGTTKSNAGMKHGTMDTESKKDPEGCDKSDDKNTVNLFGKVIETTEQDQEPSEATVGNGEVTLTYATGTKEESAGVQDNLFLEKAMQLAKRHANALFDYAVTGDVKMLLAVQRHLTAVQDENGDSVLHLAIIHLHSQLVRDLLEVTSGLISDDIINMRNDLYQTPLHLAVITKQEDVVEDLLRAGADLSLLDRLGNSVLHLAAKEGHDKVLSILLKHKKAALLLDHPNGDGLNAIHLAMMSNSLPCLLLLVAAGADVNAQEQKSGRTALHLAVEHDNISLAGCLLLEGDAHVDSTTYDGTTPLHIAAGRGSTRLAALLKAAGADPLVENFEPLYDLDDSWENAGEDEGVVPGTTPLDMATSWQVFDILNGKPYEPEFTSDDLLAQGDMKQLAEDVKLQLYKLLEIPDPDKNWATLAQKLGLGILNNAFRLSPAPSKTLMDNYEVSGGTVRELVEALRQMGYTEAIEVIQAASSPVKTTSQAHSLPLSPASTRQQIDELRDSDSVCDSGVETSFRKLSFTESLTSGASLLTLNKMPHDYGQEGPLEGKIXXMPRCPAGAMDEGPVDLRTRPKAAGLPGAALPLRKRPLRAPSPEPAAPRGAAGLVVPLDPLRGGCDLPAVPGPPHGLARPEALYYPGALLPLYPTRAMGSPFPLVNLPTPLYPMMCPMEHPLSADIAMATRADEDGDTPLHIAVVQGNLPAVHRLVNLFQQGGRELDIYNNLRQTPLHLAVITTLPSVVRLLVTAGASPMALDRHGQTAAHLACEHRSPTCLRALLDSAAPGTLDLEARNYDGLTALHVAVNTECQETVQLLLERGADIDAVDIKSGRSPLIHAVENNSLSMVQLLLQHGANVNAQMYSGSSALHSASGRGLLPLVRTLVRSGADSSLKNCHNDTPLMVARSRRVIDILRGKATRPASTSQPDPSPDRSANTSPESSSRLSSNGLLSASPSSSPSQSPPRDPPGFPMAPPNFFLPSPSPPAFLPFAGVLRGPGRPVPPSPAPGGSXXMAEDDPYLGRPEQMFHLDPSLTHTIFNPEVFQPQMALPTDGPYLQILEQPKQRGFRFRYVCEGPSHGGLPGASSEKNKKSYPQVKICNYVGPAKVIVQLVTNGKNIHLHAHSLVGKHCEDGICTVTAGPKDMVVGFANLGILHVTKKKVFETLEARMTEACIRGYNPGLLVHPDLAYLQAEGGGDRQLGDREKELIRQAALQQTKEMDLSVVRLMFTAFLPDSTGSFTRRLEPVVSDAIYDSKAPNASNLKIVRMDRTAGCVTGGEEIYLLCDKVQKDDIQIRFYEEEENGGVWEGFGDFSPTDVHRQFAIVFKTPKYKDINITKPASVFVQLRRKSDLETSEPKPFLYYPEIKDKEEVQRKRQKLMPNFSDSFGGGSGAGAGGGGMFGSGGGGGGTGSTGPGYSFPHYGFPTYGGITFHPGTTKSNAGMKHGTMDTESKKDPEGCDKSDDKNTVNLFGKVIETTEQDQEPSEATVGNGEVTLTYATGTKEESAGVQDNLFLEKAMQLAKRHANALFDYAVTGDVKMLLAVQRHLTAVQDENGDSVLHLAIIHLHSQLVRDLLEVTSGLISDDIINMRNDLYQTPLHLAVITKQEDVVEDLLRAGADLSLLDRLGNSVLHLAAKEGHDKVLSILLKHKKAALLLDHPNGDGLNAIHLAMMSNSLPCLLLLVAAGADVNAQEQKSGRTALHLAVEHDNISLAGCLLLEGDAHVDSTTYDGTTPLHIAAGRGSTRLAALLKAAGADPLVENFEPLYDLDDSWENAGEDEGVVPGTTPLDMATSWQVFDILNGKPYEPEFTSDDLLAQGDMKQLAEDVKLQLYKLLEIPDPDKNWATLAQKLGLGILNNAFRLSPAPSKTLMDNYEVSGGTVRELVEALRQMGYTEAIEVIQAASSPVKTTSQAHSLPLSPASTRQQIDELRDSDSVCDSGVETSFRKLSFTESLTSGASLLTLNKMPHDYGQEGPLEGKIX	MAEDDPYLGRPEQMFHLILL	NA	NA	NA	NA	NA	NA	NA
##	HLA-DPA10103-DPB10201	6	0_atggcagaag	EQMFHLILL	PYLGRPEQMFHLILL	285.69	14	0		1	1	chr4;chr4;chr19;chr4	NM_001165412;NM_003998;NM_005178;NM_001319226	NFKB1;NFKB1;BCL3;NFKB1	103422515;103422515;45251964;103423090	103538459;103538459;45263301;103538459	0_1	27	1	MAEDDPYLGRPEQMFHLDPSLTHTIFNPEVFQPQMALPTDGPYLQILEQPKQRGFRFRYVCEGPSHGGLPGASSEKNKKSYPQVKICNYVGPAKVIVQLVTNGKNIHLHAHSLVGKHCEDGICTVTAGPKDMVVGFANLGILHVTKKKVFETLEARMTEACIRGYNPGLLVHPDLAYLQAEGGGDRQLGDREKELIRQAALQQTKEMDLSVVRLMFTAFLPDSTGSFTRRLEPVVSDAIYDSKAPNASNLKIVRMDRTAGCVTGGEEIYLLCDKVQKDDIQIRFYEEEENGGVWEGFGDFSPTDVHRQFAIVFKTPKYKDINITKPASVFVQLRRKSDLETSEPKPFLYYPEIKDKEEVQRKRQKLMPNFSDSFGGGSGAGAGGGGMFGSGGGGGGTGSTGPGYSFPHYGFPTYGGITFHPGTTKSNAGMKHGTMDTESKKDPEGCDKSDDKNTVNLFGKVIETTEQDQEPSEATVGNGEVTLTYATGTKEESAGVQDNLFLEKAMQLAKRHANALFDYAVTGDVKMLLAVQRHLTAVQDENGDSVLHLAIIHLHSQLVRDLLEVTSGLISDDIINMRNDLYQTPLHLAVITKQEDVVEDLLRAGADLSLLDRLGNSVLHLAAKEGHDKVLSILLKHKKAALLLDHPNGDGLNAIHLAMMSNSLPCLLLLVAAGADVNAQEQKSGRTALHLAVEHDNISLAGCLLLEGDAHVDSTTYDGTTPLHIAAGRGSTRLAALLKAAGADPLVENFEPLYDLDDSWENAGEDEGVVPGTTPLDMATSWQVFDILNGKPYEPEFTSDDLLAQGDMKQLAEDVKLQLYKLLEIPDPDKNWATLAQKLGLGILNNAFRLSPAPSKTLMDNYEVSGGTVRELVEALRQMGYTEAIEVIQAASSPVKTTSQAHSLPLSPASTRQQIDELRDSDSVCDSGVETSFRKLSFTESLTSGASLLTLNKMPHDYGQEGPLEGKIXXMAEDDPYLGRPEQMFHLDPSLTHTIFNPEVFQPQMALPTADGPYLQILEQPKQRGFRFRYVCEGPSHGGLPGASSEKNKKSYPQVKICNYVGPAKVIVQLVTNGKNIHLHAHSLVGKHCEDGICTVTAGPKDMVVGFANLGILHVTKKKVFETLEARMTEACIRGYNPGLLVHPDLAYLQAEGGGDRQLGDREKELIRQAALQQTKEMDLSVVRLMFTAFLPDSTGSFTRRLEPVVSDAIYDSKAPNASNLKIVRMDRTAGCVTGGEEIYLLCDKVQKDDIQIRFYEEEENGGVWEGFGDFSPTDVHRQFAIVFKTPKYKDINITKPASVFVQLRRKSDLETSEPKPFLYYPEIKDKEEVQRKRQKLMPNFSDSFGGGSGAGAGGGGMFGSGGGGGGTGSTGPGYSFPHYGFPTYGGITFHPGTTKSNAGMKHGTMDTESKKDPEGCDKSDDKNTVNLFGKVIETTEQDQEPSEATVGNGEVTLTYATGTKEESAGVQDNLFLEKAMQLAKRHANALFDYAVTGDVKMLLAVQRHLTAVQDENGDSVLHLAIIHLHSQLVRDLLEVTSGLISDDIINMRNDLYQTPLHLAVITKQEDVVEDLLRAGADLSLLDRLGNSVLHLAAKEGHDKVLSILLKHKKAALLLDHPNGDGLNAIHLAMMSNSLPCLLLLVAAGADVNAQEQKSGRTALHLAVEHDNISLAGCLLLEGDAHVDSTTYDGTTPLHIAAGRGSTRLAALLKAAGADPLVENFEPLYDLDDSWENAGEDEGVVPGTTPLDMATSWQVFDILNGKPYEPEFTSDDLLAQGDMKQLAEDVKLQLYKLLEIPDPDKNWATLAQKLGLGILNNAFRLSPAPSKTLMDNYEVSGGTVRELVEALRQMGYTEAIEVIQAASSPVKTTSQAHSLPLSPASTRQQIDELRDSDSVCDSGVETSFRKLSFTESLTSGASLLTLNKMPHDYGQEGPLEGKIXXMPRCPAGAMDEGPVDLRTRPKAAGLPGAALPLRKRPLRAPSPEPAAPRGAAGLVVPLDPLRGGCDLPAVPGPPHGLARPEALYYPGALLPLYPTRAMGSPFPLVNLPTPLYPMMCPMEHPLSADIAMATRADEDGDTPLHIAVVQGNLPAVHRLVNLFQQGGRELDIYNNLRQTPLHLAVITTLPSVVRLLVTAGASPMALDRHGQTAAHLACEHRSPTCLRALLDSAAPGTLDLEARNYDGLTALHVAVNTECQETVQLLLERGADIDAVDIKSGRSPLIHAVENNSLSMVQLLLQHGANVNAQMYSGSSALHSASGRGLLPLVRTLVRSGADSSLKNCHNDTPLMVARSRRVIDILRGKATRPASTSQPDPSPDRSANTSPESSSRLSSNGLLSASPSSSPSQSPPRDPPGFPMAPPNFFLPSPSPPAFLPFAGVLRGPGRPVPPSPAPGGSXXMAEDDPYLGRPEQMFHLDPSLTHTIFNPEVFQPQMALPTDGPYLQILEQPKQRGFRFRYVCEGPSHGGLPGASSEKNKKSYPQVKICNYVGPAKVIVQLVTNGKNIHLHAHSLVGKHCEDGICTVTAGPKDMVVGFANLGILHVTKKKVFETLEARMTEACIRGYNPGLLVHPDLAYLQAEGGGDRQLGDREKELIRQAALQQTKEMDLSVVRLMFTAFLPDSTGSFTRRLEPVVSDAIYDSKAPNASNLKIVRMDRTAGCVTGGEEIYLLCDKVQKDDIQIRFYEEEENGGVWEGFGDFSPTDVHRQFAIVFKTPKYKDINITKPASVFVQLRRKSDLETSEPKPFLYYPEIKDKEEVQRKRQKLMPNFSDSFGGGSGAGAGGGGMFGSGGGGGGTGSTGPGYSFPHYGFPTYGGITFHPGTTKSNAGMKHGTMDTESKKDPEGCDKSDDKNTVNLFGKVIETTEQDQEPSEATVGNGEVTLTYATGTKEESAGVQDNLFLEKAMQLAKRHANALFDYAVTGDVKMLLAVQRHLTAVQDENGDSVLHLAIIHLHSQLVRDLLEVTSGLISDDIINMRNDLYQTPLHLAVITKQEDVVEDLLRAGADLSLLDRLGNSVLHLAAKEGHDKVLSILLKHKKAALLLDHPNGDGLNAIHLAMMSNSLPCLLLLVAAGADVNAQEQKSGRTALHLAVEHDNISLAGCLLLEGDAHVDSTTYDGTTPLHIAAGRGSTRLAALLKAAGADPLVENFEPLYDLDDSWENAGEDEGVVPGTTPLDMATSWQVFDILNGKPYEPEFTSDDLLAQGDMKQLAEDVKLQLYKLLEIPDPDKNWATLAQKLGLGILNNAFRLSPAPSKTLMDNYEVSGGTVRELVEALRQMGYTEAIEVIQAASSPVKTTSQAHSLPLSPASTRQQIDELRDSDSVCDSGVETSFRKLSFTESLTSGASLLTLNKMPHDYGQEGPLEGKIX	MAEDDPYLGRPEQMFHLILL	NA	NA	NA	NA	NA	NA	NA
##	DRB1_1302	3	0_atggcagaag	YLGRPEQMF	DDPYLGRPEQMFHLI	5199.08	90	0		1	1	chr4;chr4;chr19;chr4	NM_001165412;NM_003998;NM_005178;NM_001319226	NFKB1;NFKB1;BCL3;NFKB1	103422515;103422515;45251964;103423090	103538459;103538459;45263301;103538459	0_1	27	1	MAEDDPYLGRPEQMFHLDPSLTHTIFNPEVFQPQMALPTDGPYLQILEQPKQRGFRFRYVCEGPSHGGLPGASSEKNKKSYPQVKICNYVGPAKVIVQLVTNGKNIHLHAHSLVGKHCEDGICTVTAGPKDMVVGFANLGILHVTKKKVFETLEARMTEACIRGYNPGLLVHPDLAYLQAEGGGDRQLGDREKELIRQAALQQTKEMDLSVVRLMFTAFLPDSTGSFTRRLEPVVSDAIYDSKAPNASNLKIVRMDRTAGCVTGGEEIYLLCDKVQKDDIQIRFYEEEENGGVWEGFGDFSPTDVHRQFAIVFKTPKYKDINITKPASVFVQLRRKSDLETSEPKPFLYYPEIKDKEEVQRKRQKLMPNFSDSFGGGSGAGAGGGGMFGSGGGGGGTGSTGPGYSFPHYGFPTYGGITFHPGTTKSNAGMKHGTMDTESKKDPEGCDKSDDKNTVNLFGKVIETTEQDQEPSEATVGNGEVTLTYATGTKEESAGVQDNLFLEKAMQLAKRHANALFDYAVTGDVKMLLAVQRHLTAVQDENGDSVLHLAIIHLHSQLVRDLLEVTSGLISDDIINMRNDLYQTPLHLAVITKQEDVVEDLLRAGADLSLLDRLGNSVLHLAAKEGHDKVLSILLKHKKAALLLDHPNGDGLNAIHLAMMSNSLPCLLLLVAAGADVNAQEQKSGRTALHLAVEHDNISLAGCLLLEGDAHVDSTTYDGTTPLHIAAGRGSTRLAALLKAAGADPLVENFEPLYDLDDSWENAGEDEGVVPGTTPLDMATSWQVFDILNGKPYEPEFTSDDLLAQGDMKQLAEDVKLQLYKLLEIPDPDKNWATLAQKLGLGILNNAFRLSPAPSKTLMDNYEVSGGTVRELVEALRQMGYTEAIEVIQAASSPVKTTSQAHSLPLSPASTRQQIDELRDSDSVCDSGVETSFRKLSFTESLTSGASLLTLNKMPHDYGQEGPLEGKIXXMAEDDPYLGRPEQMFHLDPSLTHTIFNPEVFQPQMALPTADGPYLQILEQPKQRGFRFRYVCEGPSHGGLPGASSEKNKKSYPQVKICNYVGPAKVIVQLVTNGKNIHLHAHSLVGKHCEDGICTVTAGPKDMVVGFANLGILHVTKKKVFETLEARMTEACIRGYNPGLLVHPDLAYLQAEGGGDRQLGDREKELIRQAALQQTKEMDLSVVRLMFTAFLPDSTGSFTRRLEPVVSDAIYDSKAPNASNLKIVRMDRTAGCVTGGEEIYLLCDKVQKDDIQIRFYEEEENGGVWEGFGDFSPTDVHRQFAIVFKTPKYKDINITKPASVFVQLRRKSDLETSEPKPFLYYPEIKDKEEVQRKRQKLMPNFSDSFGGGSGAGAGGGGMFGSGGGGGGTGSTGPGYSFPHYGFPTYGGITFHPGTTKSNAGMKHGTMDTESKKDPEGCDKSDDKNTVNLFGKVIETTEQDQEPSEATVGNGEVTLTYATGTKEESAGVQDNLFLEKAMQLAKRHANALFDYAVTGDVKMLLAVQRHLTAVQDENGDSVLHLAIIHLHSQLVRDLLEVTSGLISDDIINMRNDLYQTPLHLAVITKQEDVVEDLLRAGADLSLLDRLGNSVLHLAAKEGHDKVLSILLKHKKAALLLDHPNGDGLNAIHLAMMSNSLPCLLLLVAAGADVNAQEQKSGRTALHLAVEHDNISLAGCLLLEGDAHVDSTTYDGTTPLHIAAGRGSTRLAALLKAAGADPLVENFEPLYDLDDSWENAGEDEGVVPGTTPLDMATSWQVFDILNGKPYEPEFTSDDLLAQGDMKQLAEDVKLQLYKLLEIPDPDKNWATLAQKLGLGILNNAFRLSPAPSKTLMDNYEVSGGTVRELVEALRQMGYTEAIEVIQAASSPVKTTSQAHSLPLSPASTRQQIDELRDSDSVCDSGVETSFRKLSFTESLTSGASLLTLNKMPHDYGQEGPLEGKIXXMPRCPAGAMDEGPVDLRTRPKAAGLPGAALPLRKRPLRAPSPEPAAPRGAAGLVVPLDPLRGGCDLPAVPGPPHGLARPEALYYPGALLPLYPTRAMGSPFPLVNLPTPLYPMMCPMEHPLSADIAMATRADEDGDTPLHIAVVQGNLPAVHRLVNLFQQGGRELDIYNNLRQTPLHLAVITTLPSVVRLLVTAGASPMALDRHGQTAAHLACEHRSPTCLRALLDSAAPGTLDLEARNYDGLTALHVAVNTECQETVQLLLERGADIDAVDIKSGRSPLIHAVENNSLSMVQLLLQHGANVNAQMYSGSSALHSASGRGLLPLVRTLVRSGADSSLKNCHNDTPLMVARSRRVIDILRGKATRPASTSQPDPSPDRSANTSPESSSRLSSNGLLSASPSSSPSQSPPRDPPGFPMAPPNFFLPSPSPPAFLPFAGVLRGPGRPVPPSPAPGGSXXMAEDDPYLGRPEQMFHLDPSLTHTIFNPEVFQPQMALPTDGPYLQILEQPKQRGFRFRYVCEGPSHGGLPGASSEKNKKSYPQVKICNYVGPAKVIVQLVTNGKNIHLHAHSLVGKHCEDGICTVTAGPKDMVVGFANLGILHVTKKKVFETLEARMTEACIRGYNPGLLVHPDLAYLQAEGGGDRQLGDREKELIRQAALQQTKEMDLSVVRLMFTAFLPDSTGSFTRRLEPVVSDAIYDSKAPNASNLKIVRMDRTAGCVTGGEEIYLLCDKVQKDDIQIRFYEEEENGGVWEGFGDFSPTDVHRQFAIVFKTPKYKDINITKPASVFVQLRRKSDLETSEPKPFLYYPEIKDKEEVQRKRQKLMPNFSDSFGGGSGAGAGGGGMFGSGGGGGGTGSTGPGYSFPHYGFPTYGGITFHPGTTKSNAGMKHGTMDTESKKDPEGCDKSDDKNTVNLFGKVIETTEQDQEPSEATVGNGEVTLTYATGTKEESAGVQDNLFLEKAMQLAKRHANALFDYAVTGDVKMLLAVQRHLTAVQDENGDSVLHLAIIHLHSQLVRDLLEVTSGLISDDIINMRNDLYQTPLHLAVITKQEDVVEDLLRAGADLSLLDRLGNSVLHLAAKEGHDKVLSILLKHKKAALLLDHPNGDGLNAIHLAMMSNSLPCLLLLVAAGADVNAQEQKSGRTALHLAVEHDNISLAGCLLLEGDAHVDSTTYDGTTPLHIAAGRGSTRLAALLKAAGADPLVENFEPLYDLDDSWENAGEDEGVVPGTTPLDMATSWQVFDILNGKPYEPEFTSDDLLAQGDMKQLAEDVKLQLYKLLEIPDPDKNWATLAQKLGLGILNNAFRLSPAPSKTLMDNYEVSGGTVRELVEALRQMGYTEAIEVIQAASSPVKTTSQAHSLPLSPASTRQQIDELRDSDSVCDSGVETSFRKLSFTESLTSGASLLTLNKMPHDYGQEGPLEGKIX	MAEDDPYLGRPEQMFHLILL	NA	NA	NA	NA	NA	NA	NA
##	DRB1_1302	4	0_atggcagaag	GRPEQMFHL	DPYLGRPEQMFHLIL	3820.74	80	0		1	1	chr4;chr4;chr19;chr4	NM_001165412;NM_003998;NM_005178;NM_001319226	NFKB1;NFKB1;BCL3;NFKB1	103422515;103422515;45251964;103423090	103538459;103538459;45263301;103538459	0_1	27	1	MAEDDPYLGRPEQMFHLDPSLTHTIFNPEVFQPQMALPTDGPYLQILEQPKQRGFRFRYVCEGPSHGGLPGASSEKNKKSYPQVKICNYVGPAKVIVQLVTNGKNIHLHAHSLVGKHCEDGICTVTAGPKDMVVGFANLGILHVTKKKVFETLEARMTEACIRGYNPGLLVHPDLAYLQAEGGGDRQLGDREKELIRQAALQQTKEMDLSVVRLMFTAFLPDSTGSFTRRLEPVVSDAIYDSKAPNASNLKIVRMDRTAGCVTGGEEIYLLCDKVQKDDIQIRFYEEEENGGVWEGFGDFSPTDVHRQFAIVFKTPKYKDINITKPASVFVQLRRKSDLETSEPKPFLYYPEIKDKEEVQRKRQKLMPNFSDSFGGGSGAGAGGGGMFGSGGGGGGTGSTGPGYSFPHYGFPTYGGITFHPGTTKSNAGMKHGTMDTESKKDPEGCDKSDDKNTVNLFGKVIETTEQDQEPSEATVGNGEVTLTYATGTKEESAGVQDNLFLEKAMQLAKRHANALFDYAVTGDVKMLLAVQRHLTAVQDENGDSVLHLAIIHLHSQLVRDLLEVTSGLISDDIINMRNDLYQTPLHLAVITKQEDVVEDLLRAGADLSLLDRLGNSVLHLAAKEGHDKVLSILLKHKKAALLLDHPNGDGLNAIHLAMMSNSLPCLLLLVAAGADVNAQEQKSGRTALHLAVEHDNISLAGCLLLEGDAHVDSTTYDGTTPLHIAAGRGSTRLAALLKAAGADPLVENFEPLYDLDDSWENAGEDEGVVPGTTPLDMATSWQVFDILNGKPYEPEFTSDDLLAQGDMKQLAEDVKLQLYKLLEIPDPDKNWATLAQKLGLGILNNAFRLSPAPSKTLMDNYEVSGGTVRELVEALRQMGYTEAIEVIQAASSPVKTTSQAHSLPLSPASTRQQIDELRDSDSVCDSGVETSFRKLSFTESLTSGASLLTLNKMPHDYGQEGPLEGKIXXMAEDDPYLGRPEQMFHLDPSLTHTIFNPEVFQPQMALPTADGPYLQILEQPKQRGFRFRYVCEGPSHGGLPGASSEKNKKSYPQVKICNYVGPAKVIVQLVTNGKNIHLHAHSLVGKHCEDGICTVTAGPKDMVVGFANLGILHVTKKKVFETLEARMTEACIRGYNPGLLVHPDLAYLQAEGGGDRQLGDREKELIRQAALQQTKEMDLSVVRLMFTAFLPDSTGSFTRRLEPVVSDAIYDSKAPNASNLKIVRMDRTAGCVTGGEEIYLLCDKVQKDDIQIRFYEEEENGGVWEGFGDFSPTDVHRQFAIVFKTPKYKDINITKPASVFVQLRRKSDLETSEPKPFLYYPEIKDKEEVQRKRQKLMPNFSDSFGGGSGAGAGGGGMFGSGGGGGGTGSTGPGYSFPHYGFPTYGGITFHPGTTKSNAGMKHGTMDTESKKDPEGCDKSDDKNTVNLFGKVIETTEQDQEPSEATVGNGEVTLTYATGTKEESAGVQDNLFLEKAMQLAKRHANALFDYAVTGDVKMLLAVQRHLTAVQDENGDSVLHLAIIHLHSQLVRDLLEVTSGLISDDIINMRNDLYQTPLHLAVITKQEDVVEDLLRAGADLSLLDRLGNSVLHLAAKEGHDKVLSILLKHKKAALLLDHPNGDGLNAIHLAMMSNSLPCLLLLVAAGADVNAQEQKSGRTALHLAVEHDNISLAGCLLLEGDAHVDSTTYDGTTPLHIAAGRGSTRLAALLKAAGADPLVENFEPLYDLDDSWENAGEDEGVVPGTTPLDMATSWQVFDILNGKPYEPEFTSDDLLAQGDMKQLAEDVKLQLYKLLEIPDPDKNWATLAQKLGLGILNNAFRLSPAPSKTLMDNYEVSGGTVRELVEALRQMGYTEAIEVIQAASSPVKTTSQAHSLPLSPASTRQQIDELRDSDSVCDSGVETSFRKLSFTESLTSGASLLTLNKMPHDYGQEGPLEGKIXXMPRCPAGAMDEGPVDLRTRPKAAGLPGAALPLRKRPLRAPSPEPAAPRGAAGLVVPLDPLRGGCDLPAVPGPPHGLARPEALYYPGALLPLYPTRAMGSPFPLVNLPTPLYPMMCPMEHPLSADIAMATRADEDGDTPLHIAVVQGNLPAVHRLVNLFQQGGRELDIYNNLRQTPLHLAVITTLPSVVRLLVTAGASPMALDRHGQTAAHLACEHRSPTCLRALLDSAAPGTLDLEARNYDGLTALHVAVNTECQETVQLLLERGADIDAVDIKSGRSPLIHAVENNSLSMVQLLLQHGANVNAQMYSGSSALHSASGRGLLPLVRTLVRSGADSSLKNCHNDTPLMVARSRRVIDILRGKATRPASTSQPDPSPDRSANTSPESSSRLSSNGLLSASPSSSPSQSPPRDPPGFPMAPPNFFLPSPSPPAFLPFAGVLRGPGRPVPPSPAPGGSXXMAEDDPYLGRPEQMFHLDPSLTHTIFNPEVFQPQMALPTDGPYLQILEQPKQRGFRFRYVCEGPSHGGLPGASSEKNKKSYPQVKICNYVGPAKVIVQLVTNGKNIHLHAHSLVGKHCEDGICTVTAGPKDMVVGFANLGILHVTKKKVFETLEARMTEACIRGYNPGLLVHPDLAYLQAEGGGDRQLGDREKELIRQAALQQTKEMDLSVVRLMFTAFLPDSTGSFTRRLEPVVSDAIYDSKAPNASNLKIVRMDRTAGCVTGGEEIYLLCDKVQKDDIQIRFYEEEENGGVWEGFGDFSPTDVHRQFAIVFKTPKYKDINITKPASVFVQLRRKSDLETSEPKPFLYYPEIKDKEEVQRKRQKLMPNFSDSFGGGSGAGAGGGGMFGSGGGGGGTGSTGPGYSFPHYGFPTYGGITFHPGTTKSNAGMKHGTMDTESKKDPEGCDKSDDKNTVNLFGKVIETTEQDQEPSEATVGNGEVTLTYATGTKEESAGVQDNLFLEKAMQLAKRHANALFDYAVTGDVKMLLAVQRHLTAVQDENGDSVLHLAIIHLHSQLVRDLLEVTSGLISDDIINMRNDLYQTPLHLAVITKQEDVVEDLLRAGADLSLLDRLGNSVLHLAAKEGHDKVLSILLKHKKAALLLDHPNGDGLNAIHLAMMSNSLPCLLLLVAAGADVNAQEQKSGRTALHLAVEHDNISLAGCLLLEGDAHVDSTTYDGTTPLHIAAGRGSTRLAALLKAAGADPLVENFEPLYDLDDSWENAGEDEGVVPGTTPLDMATSWQVFDILNGKPYEPEFTSDDLLAQGDMKQLAEDVKLQLYKLLEIPDPDKNWATLAQKLGLGILNNAFRLSPAPSKTLMDNYEVSGGTVRELVEALRQMGYTEAIEVIQAASSPVKTTSQAHSLPLSPASTRQQIDELRDSDSVCDSGVETSFRKLSFTESLTSGASLLTLNKMPHDYGQEGPLEGKIX	MAEDDPYLGRPEQMFHLILL	NA	NA	NA	NA	NA	NA	NA
##	DRB1_1302	3	0_atggcagaag	GRPEQMFHL	PYLGRPEQMFHLILL	3239.61	80	0		1	1	chr4;chr4;chr19;chr4	NM_001165412;NM_003998;NM_005178;NM_001319226	NFKB1;NFKB1;BCL3;NFKB1	103422515;103422515;45251964;103423090	103538459;103538459;45263301;103538459	0_1	27	1	MAEDDPYLGRPEQMFHLDPSLTHTIFNPEVFQPQMALPTDGPYLQILEQPKQRGFRFRYVCEGPSHGGLPGASSEKNKKSYPQVKICNYVGPAKVIVQLVTNGKNIHLHAHSLVGKHCEDGICTVTAGPKDMVVGFANLGILHVTKKKVFETLEARMTEACIRGYNPGLLVHPDLAYLQAEGGGDRQLGDREKELIRQAALQQTKEMDLSVVRLMFTAFLPDSTGSFTRRLEPVVSDAIYDSKAPNASNLKIVRMDRTAGCVTGGEEIYLLCDKVQKDDIQIRFYEEEENGGVWEGFGDFSPTDVHRQFAIVFKTPKYKDINITKPASVFVQLRRKSDLETSEPKPFLYYPEIKDKEEVQRKRQKLMPNFSDSFGGGSGAGAGGGGMFGSGGGGGGTGSTGPGYSFPHYGFPTYGGITFHPGTTKSNAGMKHGTMDTESKKDPEGCDKSDDKNTVNLFGKVIETTEQDQEPSEATVGNGEVTLTYATGTKEESAGVQDNLFLEKAMQLAKRHANALFDYAVTGDVKMLLAVQRHLTAVQDENGDSVLHLAIIHLHSQLVRDLLEVTSGLISDDIINMRNDLYQTPLHLAVITKQEDVVEDLLRAGADLSLLDRLGNSVLHLAAKEGHDKVLSILLKHKKAALLLDHPNGDGLNAIHLAMMSNSLPCLLLLVAAGADVNAQEQKSGRTALHLAVEHDNISLAGCLLLEGDAHVDSTTYDGTTPLHIAAGRGSTRLAALLKAAGADPLVENFEPLYDLDDSWENAGEDEGVVPGTTPLDMATSWQVFDILNGKPYEPEFTSDDLLAQGDMKQLAEDVKLQLYKLLEIPDPDKNWATLAQKLGLGILNNAFRLSPAPSKTLMDNYEVSGGTVRELVEALRQMGYTEAIEVIQAASSPVKTTSQAHSLPLSPASTRQQIDELRDSDSVCDSGVETSFRKLSFTESLTSGASLLTLNKMPHDYGQEGPLEGKIXXMAEDDPYLGRPEQMFHLDPSLTHTIFNPEVFQPQMALPTADGPYLQILEQPKQRGFRFRYVCEGPSHGGLPGASSEKNKKSYPQVKICNYVGPAKVIVQLVTNGKNIHLHAHSLVGKHCEDGICTVTAGPKDMVVGFANLGILHVTKKKVFETLEARMTEACIRGYNPGLLVHPDLAYLQAEGGGDRQLGDREKELIRQAALQQTKEMDLSVVRLMFTAFLPDSTGSFTRRLEPVVSDAIYDSKAPNASNLKIVRMDRTAGCVTGGEEIYLLCDKVQKDDIQIRFYEEEENGGVWEGFGDFSPTDVHRQFAIVFKTPKYKDINITKPASVFVQLRRKSDLETSEPKPFLYYPEIKDKEEVQRKRQKLMPNFSDSFGGGSGAGAGGGGMFGSGGGGGGTGSTGPGYSFPHYGFPTYGGITFHPGTTKSNAGMKHGTMDTESKKDPEGCDKSDDKNTVNLFGKVIETTEQDQEPSEATVGNGEVTLTYATGTKEESAGVQDNLFLEKAMQLAKRHANALFDYAVTGDVKMLLAVQRHLTAVQDENGDSVLHLAIIHLHSQLVRDLLEVTSGLISDDIINMRNDLYQTPLHLAVITKQEDVVEDLLRAGADLSLLDRLGNSVLHLAAKEGHDKVLSILLKHKKAALLLDHPNGDGLNAIHLAMMSNSLPCLLLLVAAGADVNAQEQKSGRTALHLAVEHDNISLAGCLLLEGDAHVDSTTYDGTTPLHIAAGRGSTRLAALLKAAGADPLVENFEPLYDLDDSWENAGEDEGVVPGTTPLDMATSWQVFDILNGKPYEPEFTSDDLLAQGDMKQLAEDVKLQLYKLLEIPDPDKNWATLAQKLGLGILNNAFRLSPAPSKTLMDNYEVSGGTVRELVEALRQMGYTEAIEVIQAASSPVKTTSQAHSLPLSPASTRQQIDELRDSDSVCDSGVETSFRKLSFTESLTSGASLLTLNKMPHDYGQEGPLEGKIXXMPRCPAGAMDEGPVDLRTRPKAAGLPGAALPLRKRPLRAPSPEPAAPRGAAGLVVPLDPLRGGCDLPAVPGPPHGLARPEALYYPGALLPLYPTRAMGSPFPLVNLPTPLYPMMCPMEHPLSADIAMATRADEDGDTPLHIAVVQGNLPAVHRLVNLFQQGGRELDIYNNLRQTPLHLAVITTLPSVVRLLVTAGASPMALDRHGQTAAHLACEHRSPTCLRALLDSAAPGTLDLEARNYDGLTALHVAVNTECQETVQLLLERGADIDAVDIKSGRSPLIHAVENNSLSMVQLLLQHGANVNAQMYSGSSALHSASGRGLLPLVRTLVRSGADSSLKNCHNDTPLMVARSRRVIDILRGKATRPASTSQPDPSPDRSANTSPESSSRLSSNGLLSASPSSSPSQSPPRDPPGFPMAPPNFFLPSPSPPAFLPFAGVLRGPGRPVPPSPAPGGSXXMAEDDPYLGRPEQMFHLDPSLTHTIFNPEVFQPQMALPTDGPYLQILEQPKQRGFRFRYVCEGPSHGGLPGASSEKNKKSYPQVKICNYVGPAKVIVQLVTNGKNIHLHAHSLVGKHCEDGICTVTAGPKDMVVGFANLGILHVTKKKVFETLEARMTEACIRGYNPGLLVHPDLAYLQAEGGGDRQLGDREKELIRQAALQQTKEMDLSVVRLMFTAFLPDSTGSFTRRLEPVVSDAIYDSKAPNASNLKIVRMDRTAGCVTGGEEIYLLCDKVQKDDIQIRFYEEEENGGVWEGFGDFSPTDVHRQFAIVFKTPKYKDINITKPASVFVQLRRKSDLETSEPKPFLYYPEIKDKEEVQRKRQKLMPNFSDSFGGGSGAGAGGGGMFGSGGGGGGTGSTGPGYSFPHYGFPTYGGITFHPGTTKSNAGMKHGTMDTESKKDPEGCDKSDDKNTVNLFGKVIETTEQDQEPSEATVGNGEVTLTYATGTKEESAGVQDNLFLEKAMQLAKRHANALFDYAVTGDVKMLLAVQRHLTAVQDENGDSVLHLAIIHLHSQLVRDLLEVTSGLISDDIINMRNDLYQTPLHLAVITKQEDVVEDLLRAGADLSLLDRLGNSVLHLAAKEGHDKVLSILLKHKKAALLLDHPNGDGLNAIHLAMMSNSLPCLLLLVAAGADVNAQEQKSGRTALHLAVEHDNISLAGCLLLEGDAHVDSTTYDGTTPLHIAAGRGSTRLAALLKAAGADPLVENFEPLYDLDDSWENAGEDEGVVPGTTPLDMATSWQVFDILNGKPYEPEFTSDDLLAQGDMKQLAEDVKLQLYKLLEIPDPDKNWATLAQKLGLGILNNAFRLSPAPSKTLMDNYEVSGGTVRELVEALRQMGYTEAIEVIQAASSPVKTTSQAHSLPLSPASTRQQIDELRDSDSVCDSGVETSFRKLSFTESLTSGASLLTLNKMPHDYGQEGPLEGKIX	MAEDDPYLGRPEQMFHLILL	NA	NA	NA	NA	NA	NA	NA
```

Count the number of neoantigens for each peptide of which threshold is IC 50 of mutatnt peptides < 500. 
```
[R]
print(Export_Summary_Fragments(Result_HLA2_Seq, Mut_IC50_th = 500))
```

```
## [[1]]
##                                                     atggcagaag
## Num_Peptide_Per_Grp                                      3.000
## Num_Cond_Peptide_Per_Grp                                 3.000
## Num_Rest_Peptide_Per_Grp                                 1.000
## Num_Rest_Peptide_Per_Grp / Num_Cond_Peptide_Per_Grp      0.333
## -logP                                                    0.460
```

```
## [[2]]
##                                                        
## Num_Peptide_Per_NM                                3.000
## Num_Cond_Peptide_Per_NM                           3.000
## Num_Rest_Peptide_Per_NM                           1.000
## Num_Rest_Peptide_Per_NM / Num_Cond_Peptide_Per_NM 0.333
## -logP                                             0.460
```

```
## [[3]]
##                                                     MAEDDPYLGRPEQMFHLILL-0_atggcagaag-0_1
## Num_Peptide_Per_Pep                                                                 3.000
## Num_Cond_Peptide_Per_Pep                                                            3.000
## Num_Rest_Peptide_Per_Pep                                                            1.000
## Num_Rest_Peptide_Per_Pep / Num_Cond_Peptide_Per_Pep                                 0.333
## -logP                                                                               0.460
```