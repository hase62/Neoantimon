## Updated on 18, Apr. 2018. 
##
## 1. Preparation
**-Download and Set netMHCpan4.0 (Required)**

1. Download netMHCpan4.0 from http://www.cbs.dtu.dk/cgi-bin/nph-sw_request?netMHCpan. 

2. Download a script from https://github.com/hase62/Neoantimon/raw/master/lib/setNetMHCpan4.0.sh and run it as "./setNetMHCpan4.0.sh". 

**-Download and Set netMHCpan3.0 (If you require the old version)**

1. Download netMHCpan3.0 from http://www.cbs.dtu.dk/cgi-bin/sw_request?netMHCpan+3.0. 

2. Download a script from https://github.com/hase62/Neoantimon/raw/master/lib/setNetMHCpan3.0.sh and run it as "./setNetMHCpan3.0.sh". 

**-Download and Set netMHCIIpan (Required)**

1. Download netMHCIIpan 3.1 from http://www.cbs.dtu.dk/cgi-bin/nph-sw_request?netMHCIIpan. 

2. Download a script from https://github.com/hase62/Neoantimon/raw/master/lib/setNetMHCIIpan3.1.sh and run it as "./setNetMHCIIpan3.1.sh". 

**-Download refMrna Files (Required)**

**You have to get your corresponding version from GRCh38, hg38, GRCh37 or hg19.**

GRCh38/hg38: Download refMrna Files from "http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/refMrna.fa.gz". 
Otherwise, run the following codes or use "InstallRefMrnaFile(url = "http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/refMrna.fa.gz")" after installing Neoantimon. 
```
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/refMrna.fa.gz
gunzip refMrna.fa.gz
```

GRCh37/hg19: Download refMrna Files from "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/refMrna.fa.gz". 
Otherwise, run the following codes or use "InstallRefMrnaFile(url = "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/refMrna.fa.gz")" after installing Neoantimon. 
```
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/refMrna.fa.gz
gunzip refMrna.fa.gz
```

**-Download refFlat Files (Required)**

**You have to get your corresponding version from GRCh38, hg38, GRCh37 or hg19.**

GRCh38/hg38: Download refFlat Files from "http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/refFlat.txt.gz". 
Otherwise, run the following codes or use "InstallRefFlat(url = "http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/refFlat.txt.gz")" after installing Neoantimon. 
```
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/refFlat.txt.gz
gunzip refFlat.txt.gz
```

GRCh37/hg19: Download refFlat Files from "http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/refFlat.txt.gz". 
Otherwise, run the following codes or use "InstallRefFlat(url = "http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/refFlat.txt.gz")" after installing Neoantimon. 
```
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/refFlat.txt.gz
gunzip refFlat.txt.gz
```

**-Install Samtools 0_x_x**

**Required for SV fusions. Not Required for Snv/Indel, but please download if you want to calculate Allele Specific RNA Expression based on RNA bam.**

You can install samtools_0_x_x version from https://sourceforge.net/projects/samtools/files/samtools/0.1.19/samtools-0.1.19.tar.bz2. 
Otherwise, run the following codes or use "InstallSamtools()" after installing Neoantimon. 
```
wget https://sourceforge.net/projects/samtools/files/samtools/0.1.19/samtools-0.1.19.tar.bz2
tar jxf samtools-0.1.19.tar.bz2
cd samtools-0.1.19
make
cd ..
```

**-Download human refSeq**

**Required for SV fusions. Not Required for Snv/Indel, but please download if you want to calculate Allele Specific RNA Expression using RNA bam.**

**You have to get your corresponding version from GRCh38, hg38, GRCh37 or hg19.**

GRCh38: Download human refSeq by
```
wget ftp://ftp.ensembl.org/pub/release-87/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.toplevel.fa.gz
mv Homo_sapiens.GRCh38.dna.toplevel.fa.gz GRCh38.fa.gz
gunzip GRCh38.fa.gz
samtools faidx GRCh38.fa
```

hg38: Download human refSeq by:
```
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
gunzip hg38.fa.gz
samtools faidx hg38.fa
```

GRCh37: Download human refSeq:
```
wget ftp://ftp.ensembl.org/pub/release-75//fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.75.dna.toplevel.fa.gz
mv Homo_sapiens.GRCh37.75.dna.toplevel.fa.gz GRCh37.fa.gz
gunzip GRCh37.fa.gz
samtools faidx GRCh37.fa
```

**-Download SampleFiles (Not Required)**

You can get these from https://github.com/hase62/Neoantimon/raw/master/lib/data.zip. 
Otherwise, run the following codes or use "InstallSampleFiles()" after installing Neoantimon. 
```
wget https://github.com/hase62/Neoantimon/raw/master/lib/data.zip
unzip data.zip
```

## 2. Use on R
```
install.packages("devtools");
library(devtools);
install_github('hase62/Neoantimon');
library(Neoantimon);
```

## 3. Data Format
**-HLA Table**

**1. A HLA Class I table file must be according to the following format. **
```r
library(Neoantimon)
data(“sample_hla_table_c1”)
print(sample_hla_table_c1, row.names = FALSE)
```

```
##     Name      A1      A2      B1      B2      C1      C2
##   sample A*02:01 A*32:01 B*15:17 B*51:01 C*07:01 C*15:02
##  sample2 A*02:01 A*32:01 B*15:17 B*51:01 C*07:01 C*15:02
##  ...
```

**2. A HLA Class II table file must be according to the following format. **
```r
data("sample_hla_table_c2”)
print(sample_hla_table_c2, row.names = FALSE)
```

```
##	Name	DPA11	DPA12	DPB11	DPB12	DQA11	DQA12	DQB11	DQB12	DRB11	DRB12
##	sample	DPA1*01:03	DPA1*02:01	DPB1*02:01	DPB1*09:01	DQA1*01:02	DQA1*05:05	DQB1*03:01	DQB1*06:04	DRB1*11:04	DRB1*13:02
##	sample2	DPA1*01:03	DPA1*02:01	DPB1*02:01	DPB1*09:01	DQA1*01:02	DQA1*05:05	DQB1*03:01	DQB1*06:04	DRB1*11:04	DRB1*13:02
##  ...
```

**-Annotated VCF file**

**An annotated VCF file is required for Snv/Indel.**

It must include columns representing "Chromosome Number", "Mutation Start Position", "Mutation End Position", "Mutation Ref", "Mutation Alt", and "NM_ID (AAChange.refGene)".
Annotations "Chr", "Start", "End", "Ref", "Alt", "AAChange.refGene", "Depth_tumor", and "Depth_normal" are automatically detected. Otherwise, you must indicate columns for them when using Main**() functions. 
```r
data(“sample_vcf”)
print(sample_vcf, row.names = FALSE)
```

```
## Chr     Start       End Ref Alt Func.refGene Gene.refGene   GeneDetail.refGene ExonicFunc.refGene                              AAChange.refGene   cytoBand depth_tumor variantNum_tumor depth_normal
##   1  47399872  47399872   A   G       exonic      CYP4A11        nonsynonymous                SNV      CYP4A11:NM_000778:exon8:c.T1064C:p.L355P       1p33          64               28           49
##   1 116941338 116941338   T   C       exonic       ATP1A1           synonymous                SNV   ATP1A1:NM_001160234:exon16:c.T2127C:p.D709D     1p13.1         100               39          111
##   4  24556416  24556416   T   C       exonic        DHX15        nonsynonymous                SNV        DHX15:NM_001358:exon5:c.A1012G:p.T338A     4p15.2         143               47          151
##   4  70156404  70156404   -   T       exonic      UGT2B28           frameshift          insertion   UGT2B28:NM_053039:exon5:c.1186dupT:p.L395fs     4q13.2          43               15           41
##   6  75899298  75899298   T   -       exonic      COL12A1           frameshift           deletion    COL12A1:NM_004370:exon6:c.628delA:p.I210fs       6q13         122               38           73
##   9  89561162  89561162   C   T       exonic         GAS1        nonsynonymous                SNV          GAS1:NM_002048:exon1:c.G533A:p.R178H    9q21.33          20                5           26
```

**-Annotated BND format VCF file**

*An annotated BND format VCF file is required for SV fusion. *

It must include columns representing "Chromosome Number", "Mutation Start Position", "Mutation End Position", "Mutation Ref", "Mutation Alt", and "NM_ID (AAChange.refGene)" or "Gene Symbol (Gene.refGene)".
Annotations "Chr", "Start", "End", "Ref", "Alt", "Depth_tumor", and "Depth_normal" are automatically detected. Otherwise, you must indicate columns for them when using Main**() functions. 
```r
data(“sample_sv_bnd”)
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
```

**-RNAseq file (Not Required)**

**An RNAseq file is not required, but you can attach "RNA expression" information by indicating "rnaexp_file".
If you also indicate "rnabam_file", variant allele frequencies are also attached. **

```r
data(“sample_rna_exp”)
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

**-CNV file (Not Required)**

**A copynumber file is not required, but you can attach "Copy Number" information by indicating "cnv_file" and "purity".
Purity is set 1 as default value. **

```r
data(“sample_copynum”)
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

Sample files can be downloaded from https://github.com/hase62/Neoantimon/raw/master/lib/data.zip. 
```
lib/data/sample_copynum.txt
lib/data/sample_hla_table_c1.txt
lib/data/sample_hla_table_c2.txt
lib/data/sample_rna_exp.txt
lib/data/sample_vcf.txt
lib/data/sample_sv_bnd.txt
```

Prepare the following files using the above explanations. 
```
lib/netMHCpan-4.0
lib/netMHCIIpan-3.1  
lib/samtools-0.1.19
lib/refFlat.txt 
lib/refMrna.fa
lib/GRCh37.fa
```

Analize sample files ... 
```
install.packages("devtools");
library(devtools);
install_github('hase62/Neoantimon');
library(Neoantimon);
```

Calculate Neoantigens on SNVs/INDELs for HLA Class I and II. 
```
  MainSNVClass1(input_file = "lib/data/sample_vcf.txt",
                file_name_in_hla_table = "sample",
                hla_file = "lib/data/sample_hla_table_c1.txt",
                refflat_file  = "lib/refFlat.txt",
                refmrna_file = "lib/refMrna.fa",
                rnaexp_file = "lib/data/sample_rna_exp.txt",
                netMHCpan_dir = "lib/netMHCpan-3.0/netMHCpan",
                nm_id_column = 10,
                depth_tumor_column = 12,
                depth_normal_column = 14
  )

  MainSNVClass2(input_file = "lib/data/sample_vcf.txt",
                file_name_in_hla_table = "sample",
                hla_file = "lib/data/sample_hla_table_c2.txt",
                refflat_file  = "lib/refFlat.txt",
                refmrna_file = "lib/refMrna.fa",
                rnaexp_file = "lib/data/sample_rna_exp.txt",
                netMHCIIpan_dir = "lib/netMHCIIpan-3.1/netMHCIIpan",
                nm_id_column = 10,
                depth_tumor_column = 12,
                depth_normal_column = 14
  )
 
  MainINDELClass1(input_file = "lib/data/sample_vcf.txt",
                  file_name_in_hla_table = "sample",
                  hla_file = "lib/data/sample_hla_table_c1.txt",
                  refflat_file  = "lib/refFlat.txt",
                  refmrna_file = "lib/refMrna.fa",
                  rnaexp_file = "lib/data/sample_rna_exp.txt",
                  netMHCpan_dir = "lib/netMHCpan-3.0/netMHCpan",
                  nm_id_column = 10,
                  depth_tumor_column = 12,
                  depth_normal_column = 14
  )

  MainINDELClass2(input_file = "lib/data/sample_vcf.txt",
                  file_name_in_hla_table = "sample",
                  hla_file = "lib/data/sample_hla_table_c2.txt",
                  refflat_file  = "lib/refFlat.txt",
                  refmrna_file = "lib/refMrna.fa",
                  rnaexp_file = "lib/data/sample_rna_exp.txt",
                  netMHCIIpan_dir = "lib/netMHCIIpan-3.1/netMHCIIpan",
                  nm_id_column = 10,
                  depth_tumor_column = 12,
                  depth_normal_column = 14
  )
```

Merge Neoantigens on SNVs/INDELs for HLA Class I and II. 
```
  MainMergeSNVClass1(input_dir = "result.sample.NO_job_id_SNV",
                     file_prefix = "NO_job_id_SNV",
                     annotation_file = "lib/data/sample_vcf.txt.NO_job_id_SNV.peptide.txt"
  )

  MainMergeSNVClass2(input_dir = "result.sample.NO_job_id_SNV",
                     file_prefix = "NO_job_id_SNV",
                     annotation_file = "lib/data/sample_vcf.txt.NO_job_id_SNV.peptide.txt"
  )

  MainMergeINDELSVClass1(input_dir = "result.sample.NO_job_id_INDEL",
                         file_prefix = "NO_job_id_INDEL",
                         annotation_file = "lib/data/sample_vcf.txt.NO_job_id_INDEL.peptide.txt"
  )

  MainMergeINDELSVClass2(input_dir = "result.sample.NO_job_id_INDEL",
                         file_prefix = "NO_job_id_INDEL",
                         annotation_file = "lib/data/sample_vcf.txt.NO_job_id_INDEL.peptide.txt"
  )
```

Calculate Neoantigens on SV fusions for HLA Class I and II. 
```
  MainSVFUSIONClass1(input_file = "lib/data/sample_sv_bnd.txt",
                     file_name_in_hla_table = "sample",
                     hla_file = "lib/data/sample_hla_table_c1.txt",
                     refflat_file  = "lib/refFlat.txt",
                     refmrna_file = "lib/refMrna.fa",
                     rnaexp_file = "lib/data/sample_rna_exp.txt",
                     netMHCpan_dir = "lib/netMHCpan-3.0/netMHCpan",
                     refdna_file = "lib/GRCh37.fa",
                     mutation_alt_bnd_column = 5,
                     gene_symbol_column = 7,
                     mate_id_column = 8
  )

  MainSVFUSIONClass2(input_file = "lib/data/sample_sv_bnd.txt",
                     file_name_in_hla_table = "sample",
                     hla_file = "lib/data/sample_hla_table_c2.txt",
                     refflat_file  = "lib/refFlat.txt",
                     refmrna_file = "lib/refMrna.fa",
                     rnaexp_file = "lib/data/sample_rna_exp.txt",
                     netMHCIIpan_dir = "lib/netMHCIIpan-3.1/netMHCIIpan",
                     refdna_file = "lib/GRCh37.fa",
                     mutation_alt_bnd_column = 5,
                     gene_symbol_column = 7,
                     mate_id_column = 8
  )
```

Merge Neoantigens on SV fusions for HLA Class I and II. 
```
  MainMergeINDELSVClass1(input_dir = "result.sample.NO_job_id_SVFusion",
                         file_prefix = "NO_job_id_SVFusion",
                         annotation_file = "lib/data/sample_sv_bnd.txt.NO_job_id_SVFusion.peptide.txt"
  )

  MainMergeINDELSVClass2(input_dir = "result.sample.NO_job_id_SVFusion",
                         file_prefix = "NO_job_id_SVFusion",
                         annotation_file = "lib/data/sample_sv_bnd.txt.NO_job_id_SVFusion.peptide.txt"
  )
```

## 5. Result

sample_result_SNV_CLASS1_ALL
```
##	HLA	Pos	Gene	MutatedPeptide	Mut_IC50	Mut_Rank	Norm_Peptide	Norm_IC50	Norm_Rank		Gene ID	Chr	NM_ID	Change	ref	alt	Prob	Mutation Prob.	Exon Start	Exon End	Mutation Position	Depth	TumorDepth	Peptide Normal	Peptide Mutation	TotalRNA	TumorRNARatio	TumorRNA	nA	nB	Checker	MutRatio	MutRatio Min	MutRatio Max
##	HLA-A*02:01	2	0_CYP4A11	HQERCREEIHSLP	37362.9	55.00	HQERCREEIHSLL	27284.3	32.00	1	0_CYP4A11	1	NM_000778	c.T1064C	A	G	0	0	47394845	47407156	1_47399872	49	113	KHQERCREEIHSLLGDGASITWNHLDQ	KHQERCREEIHSLPGDGASITWNHLDQ	NA	NA	NA	NA	NA	NA	NA	NA	NA
##	HLA-A*02:01	3	0_CYP4A11	QERCREEIHSLPG	45831.8	90.00	QERCREEIHSLLG	44161.6	85.00	1	0_CYP4A11	1	NM_000778	c.T1064C	A	G	0	0	47394845	47407156	1_47399872	49	113	KHQERCREEIHSLLGDGASITWNHLDQ	KHQERCREEIHSLPGDGASITWNHLDQ	NA	NA	NA	NA	NA	NA	NA	NA	NA
##	HLA-A*02:01	4	0_CYP4A11	ERCREEIHSLPGD	45989.7	90.00	ERCREEIHSLLGD	44575.4	85.00	1	0_CYP4A11	1	NM_000778	c.T1064C	A	G	0	0	47394845	47407156	1_47399872	49	113	KHQERCREEIHSLLGDGASITWNHLDQ	KHQERCREEIHSLPGDGASITWNHLDQ	NA	NA	NA	NA	NA	NA	NA	NA	NA
##	HLA-A*02:01	5	0_CYP4A11	RCREEIHSLPGDG	46482.0	95.00	RCREEIHSLLGDG	44377.6	85.00	1	0_CYP4A11	1	NM_000778	c.T1064C	A	G	0	0	47394845	47407156	1_47399872	49	113	KHQERCREEIHSLLGDGASITWNHLDQ	KHQERCREEIHSLPGDGASITWNHLDQ	NA	NA	NA	NA	NA	NA	NA	NA	NA
##	HLA-A*02:01	6	0_CYP4A11	CREEIHSLPGDGA	42936.1	80.00	CREEIHSLLGDGA	41760.9	70.00	1	0_CYP4A11	1	NM_000778	c.T1064C	A	G	0	0	47394845	47407156	1_47399872	49	113	KHQERCREEIHSLLGDGASITWNHLDQ	KHQERCREEIHSLPGDGASITWNHLDQ	NA	NA	NA	NA	NA	NA	NA	NA	NA
```
