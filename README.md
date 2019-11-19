# Major Version will be released soon. 
# A manuscript for this journal is at bioarxiv (). 

## 1. Preparation
### -Install wget  (Required)
```
brew install wget
```

### -Download and Set netMHCpan4.0 (Required)

1. Download netMHCpan4.0 from http://www.cbs.dtu.dk/cgi-bin/nph-sw_request?netMHCpan and move it to a working directory. 

2. Do initial setting as followings.
```
wget --no-check-certificate https://github.com/hase62/Neoantimon/raw/master/lib/setNetMHCpan4.0.sh
chmod 750 setNetMHCpan4.0.sh
./setNetMHCpan4.0.sh 
```
### -Download and Set mhcflurry (Not Required)
1. (Recommended) Install anaconda from https://www.anaconda.com/distribution/, then run the following codes. 
```
pip install mhcflurry
mhcflurry-downloads fetch
```

2. Otherwise, install python from https://www.python.org/downloads/release/python-380/, then run the above codes.

### -Download and Set netMHCIIpan3.2 (Required)

1. Download netMHCIIpan 3.2 from http://www.cbs.dtu.dk/cgi-bin/nph-sw_request?netMHCIIpan and move it to a working directory. 

2. Do initial setting as followings.
```
wget --no-check-certificate https://github.com/hase62/Neoantimon/raw/master/lib/setNetMHCIIpan3.2.sh
chmod 750 setNetMHCIIpan3.2.sh
./setNetMHCIIpan3.2.sh
```

### -Download refMrna Files (Required)

**(You have to select your corresponding version from GRCh38, hg38, GRCh37 or hg19.)**

**GRCh38/hg38**: Run the following codes. 
```
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/refMrna.fa.gz
gunzip refMrna.fa.gz
```

**GRCh37/hg19**: Run the following codes. 
```
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/refMrna.fa.gz
gunzip refMrna.fa.gz
```

### -Download refFlat Files (Required)

**(You have to select your corresponding version from GRCh38, hg38, GRCh37 or hg19.)**

**GRCh38/hg38**: Run the following codes. 
```
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/refFlat.txt.gz
gunzip refFlat.txt.gz
```

**GRCh37/hg19**: Run the following codes. 
```
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/refFlat.txt.gz
gunzip refFlat.txt.gz
```

### -Install Samtools (Not Required)

**Required for SV fusions. Not Required for Snv/Indel, but please download if you want to calculate Allele Specific RNA Expression based on RNA bam.**
1. (Recommended) Install anaconda from https://www.anaconda.com/distribution/, then run the following codes. 
```
conda install -c bioconda samtools
conda install -c bioconda/label/cf201911 samtools
```

2. Otherwise, you can install local samtools as followings. 
```
wget https://sourceforge.net/projects/samtools/files/samtools/0.1.19/samtools-0.1.19.tar.bz2
tar jxf samtools-0.1.19.tar.bz2
```

### -Download human refSeq

**Required for SV fusions. Not Required for Snv/Indel, but please download if you want to calculate Allele Specific RNA Expression using RNA bam.**

**You have to select your corresponding version from GRCh38, hg38, GRCh37 or hg19.**

**GRCh38**: Run the following codes.
```
wget ftp://ftp.ensembl.org/pub/release-87/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.toplevel.fa.gz
mv Homo_sapiens.GRCh38.dna.toplevel.fa.gz GRCh38.fa.gz
gunzip GRCh38.fa.gz
samtools faidx GRCh38.fa
```

**hg38**: Run the following codes.
```
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
gunzip hg38.fa.gz
samtools faidx hg38.fa
```

**GRCh37/hg19**: Run the following codes.
```
wget ftp://ftp.ensembl.org/pub/release-75//fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.75.dna.toplevel.fa.gz
mv Homo_sapiens.GRCh37.75.dna.toplevel.fa.gz GRCh37.fa.gz
gunzip GRCh37.fa.gz
samtools faidx GRCh37.fa
```

### -Download SampleFiles (Not required)

Run the following codes. 
```
wget --no-check-certificate https://github.com/hase62/Neoantimon/raw/master/lib/data.zip
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
### -HLA Table

**1. A HLA Class I table file must be according to the following format.**

```r
library(Neoantimon)
data("sample_hla_table_c1")
print(sample_hla_table_c1, row.names = FALSE)
```

```
##     Name      A1      A2      B1      B2      C1      C2
##   sample A*02:01 A*32:01 B*15:17 B*51:01 C*07:01 C*15:02
##  sample2 A*02:01 A*32:01 B*15:17 B*51:01 C*07:01 C*15:02
##  ...
```

**2. A HLA Class II table file must be according to the following format.**

```r
data("sample_hla_table_c2")
print(sample_hla_table_c2, row.names = FALSE)
```

```
##	Name	DPA11	DPA12	DPB11	DPB12	DQA11	DQA12	DQB11	DQB12	DRB11	DRB12
##	sample	DPA1*01:03	DPA1*02:01	DPB1*02:01	DPB1*09:01	DQA1*01:02	DQA1*05:05	DQB1*03:01	DQB1*06:04	DRB1*11:04	DRB1*13:02
##	sample2	DPA1*01:03	DPA1*02:01	DPB1*02:01	DPB1*09:01	DQA1*01:02	DQA1*05:05	DQB1*03:01	DQB1*06:04	DRB1*11:04	DRB1*13:02
##  ...
```

### -Annotated VCF file

**An annotated VCF file is required for Snv/Indel.**

It must include columns representing "Chromosome Number", "Mutation Start Position", "Mutation End Position", "Mutation Ref", "Mutation Alt", and "NM_ID (AAChange.refGene)".
Annotations "Chr", "Start", "End", "Ref", "Alt", "AAChange.refGene", "Depth_tumor", and "Depth_normal" are automatically detected. Otherwise, you have to manually indicate columns. 
```r
data("sample_vcf")
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

### -Annotated BND format VCF file

*An annotated BND format VCF file is required for SV fusion.*

It must include columns representing "Chromosome Number", "Mutation Start Position", "Mutation End Position", "Mutation Ref", "Mutation Alt", and "NM_ID (AAChange.refGene)" or "Gene Symbol (Gene.refGene)".
Annotations "Chr", "Start", "End", "Ref", "Alt", "Depth_tumor", and "Depth_normal" are automatically detected. Otherwise, you have to manually indicate columns. 
```r
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
```

### -RNA expression file

An RNA expressoin file is not required, but you can attach "RNA expression" information by indicating "rnaexp_file" in main functions.
If you also indicate "rnabam_file", variant allele frequencies and tumor specific RNA expressions are also attached to the results. 

```r
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
They are used to calculate tumor subclonal cell population 
Purity is set 1 as default value. 

```r
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

##### Sample files can be downloaded from https://github.com/hase62/Neoantimon/raw/master/lib/data.zip. 
```
lib/data/sample_result_INDEL_CLASS1_ALL.txt
lib/data/sample_result_INDEL_CLASS2_ALL.txt
lib/data/sample_result_SeqFragment_CLASS1_ALL.txt
lib/data/sample_result_SeqFragment_CLASS2_ALL.txt
lib/data/sample_result_SNV_CLASS1_ALL.txt
lib/data/sample_result_SNV_CLASS2_ALL.txt
lib/data/sample_result_SVFusion_CLASS1_ALL.txt
lib/data/sample_result_SVFusion_CLASS2_ALL.txt
lib/data/sample_copynum.txt
lib/data/sample_hla_table_c1.txt
lib/data/sample_hla_table_c2.txt
lib/data/sample_rna_exp.txt
lib/data/sample_vcf.txt
lib/data/sample_sv_bnd.txt
```

##### Prepare the following files using the above explanations. 
```
lib/netMHCpan-4.0
lib/netMHCIIpan-3.1  
lib/samtools-0.1.19
lib/refFlat.txt 
lib/refMrna.fa
lib/GRCh37.fa
```

##### Preparation
```
install.packages("devtools");
library(devtools);
install_github('hase62/Neoantimon');
library(Neoantimon);
```

##### Calculate Neoantigens on SNVs/INDELs for HLA Class I and II. 
```
  Result_HLA1_SNV <- MainSNVClass1(input_file = "lib/data/sample_vcf.txt",
                                   file_name_in_hla_table = "sample",
                                   hla_file = "lib/data/sample_hla_table_c1.txt",
                                   refflat_file  = "lib/refFlat.txt",
                                   refmrna_file = "lib/refMrna.fa",
                                   rnaexp_file = "lib/data/sample_rna_exp.txt",
                                   netMHCpan_dir = "lib/netMHCpan-4.0/netMHCpan",
                                   depth_tumor_column = 12,
                                   depth_normal_column = 14
  )

  Result_HLA2_SNV <- MainSNVClass2(input_file = "lib/data/sample_vcf.txt",
                                   file_name_in_hla_table = "sample",
                                   hla_file = "lib/data/sample_hla_table_c2.txt",
                                   refflat_file  = "lib/refFlat.txt",
                                   refmrna_file = "lib/refMrna.fa",
                                   rnaexp_file = "lib/data/sample_rna_exp.txt",
                                   netMHCIIpan_dir = "lib/netMHCIIpan-3.1/netMHCIIpan",
                                   depth_tumor_column = 12,
                                   depth_normal_column = 14
  )

  Result_HLA1_INDEL <- MainINDELClass1(input_file = "lib/data/sample_vcf.txt",
                                       file_name_in_hla_table = "sample",
                                       hla_file = "lib/data/sample_hla_table_c1.txt",
                                       refflat_file  = "lib/refFlat.txt",
                                       refmrna_file = "lib/refMrna.fa",
                                       rnaexp_file = "lib/data/sample_rna_exp.txt",
                                       netMHCpan_dir = "lib/netMHCpan-4.0/netMHCpan",
                                       depth_tumor_column = 12,
                                       depth_normal_column = 14
  )

  Result_HLA2_INDEL <- MainINDELClass2(input_file = "lib/data/sample_vcf.txt",
                                       file_name_in_hla_table = "sample",
                                       hla_file = "lib/data/sample_hla_table_c2.txt",
                                       refflat_file  = "lib/refFlat.txt",
                                       refmrna_file = "lib/refMrna.fa",
                                       rnaexp_file = "lib/data/sample_rna_exp.txt",
                                       netMHCIIpan_dir = "lib/netMHCIIpan-3.1/netMHCIIpan",
                                       depth_tumor_column = 12,
                                       depth_normal_column = 14
  )
```

##### Calculate Neoantigens on SV fusions for HLA Class I and II. 
```
  Result_HLA1_SV <- MainSVFUSIONClass1(input_file = "lib/data/sample_sv_bnd.txt",
                                       file_name_in_hla_table = "sample",
                                       hla_file = "lib/data/sample_hla_table_c1.txt",
                                       refflat_file  = "lib/refFlat.txt",
                                       refmrna_file = "lib/refMrna.fa",
                                       rnaexp_file = "lib/data/sample_rna_exp.txt",
                                       netMHCpan_dir = "lib/netMHCpan-4.0/netMHCpan",
                                       refdna_file = "lib/GRCh37.fa",
                                       mutation_alt_bnd_column = 5,
                                       gene_symbol_column = 7,
                                       mate_id_column = 8
  )

  Result_HLA2_SV <- MainSVFUSIONClass2(input_file = "lib/data/sample_sv_bnd.txt",
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

##### Calculate Neoantigens from a fragment of RNA sequence for HLA Class I and II by comparing to the original protein. 
```
  Result_HLA1_Seq <- MainSeqFragmentClass1(input_sequence = c("atggcagaagatgatccatatttgggaaggcctgaaaaaatgtttcatttggatccttctttgactcatacaatatttaatc",
  															  "tggcagaagatgatccatatttgggaaggcctgaaaaaatgtttcatttggatccttctttgactcatacaatatttaatc"),
                                           file_name_in_hla_table = "sample",
                                           hla_file = "lib/data/sample_hla_table_c1.txt",
                                           hmdir = getwd(),
                                           job_id = "NO_job_id",
                                           refflat_file  = "lib/refFlat.txt",
                                           refmrna_file = "lib/refMrna.fa",
                                           netMHCpan_dir = "lib/netMHCpan-4.0/netMHCpan",
                                           nm_id = c("NM_003998", "NM_001165412"),
                                           reading_frame = 1
  )

  Result_HLA2_Seq <- MainSeqFragmentClass2(input_nm_id = c("NM_003998", "NM_001165412"),
                                           file_name_in_hla_table = "sample",
                                           hla_file = "lib/data/sample_hla_table_c2.txt",
                                           hmdir = getwd(),
                                           job_id = "NO_job_id",
                                           refflat_file  = "lib/refFlat.txt",
                                           refmrna_file = "lib/refMrna.fa",
                                           netMHCIIpan_dir = "lib/netMHCIIpan-3.1/netMHCIIpan",
                                           gene_symbol = c("NFKB1", "BCL3"),
                                           reading_frame = 3
  )
```

## 5. Result

**They're also included in Result_HLA*_**

sample_result_SNV_CLASS1_ALL.txt
```
##           HLA Pos    Gene Evaluated_Mutant_Peptide Mut_IC50 Mut_Rank
## 1 HLA-A*02:01   2 0_DHX15            EPERDYLEAAIRA  37952.9  63.3899
## 2 HLA-A*02:01   3 0_DHX15            PERDYLEAAIRAV  12108.2  16.5100
## 3 HLA-A*02:01   4 0_DHX15            ERDYLEAAIRAVI  12109.5  16.5114
## 4 HLA-A*02:01   5 0_DHX15            RDYLEAAIRAVIQ  32063.2  46.0009
## 5 HLA-A*02:01   6 0_DHX15            DYLEAAIRAVIQI   4819.6   9.2763
## 6 HLA-A*02:01   7 0_DHX15            YLEAAIRAVIQIH   7309.7  11.8086
##   Evaluated_Wt_Peptide Wt_IC50 Wt_Rank Chr     NM_ID   Change Ref Alt Prob
## 1        EPERDYLEAAIRT 41213.3 75.7538   4 NM_001358 c.A1012G   T   C    0
## 2        PERDYLEAAIRTV 11859.7 16.2523   4 NM_001358 c.A1012G   T   C    0
## 3        ERDYLEAAIRTVI 14831.9 19.2905   4 NM_001358 c.A1012G   T   C    0
## 4        RDYLEAAIRTVIQ 33185.2 48.8110   4 NM_001358 c.A1012G   T   C    0
## 5        DYLEAAIRTVIQI  4504.5  8.9225   4 NM_001358 c.A1012G   T   C    0
## 6        YLEAAIRTVIQIH  6924.0 11.4357   4 NM_001358 c.A1012G   T   C    0
##   Mutation_Prob. Exon_Start Exon_End Mutation_Position Total_Depth
## 1              0   24529087 24586184        4_24556416         294
## 2              0   24529087 24586184        4_24556416         294
## 3              0   24529087 24586184        4_24556416         294
## 4              0   24529087 24586184        4_24556416         294
## 5              0   24529087 24586184        4_24556416         294
## 6              0   24529087 24586184        4_24556416         294
##   Tumor_Depth                  Wt_Peptide              Mutant_Peptide
## 1         143 PEPERDYLEAAIRTVIQIHMCEEEEGD PEPERDYLEAAIRAVIQIHMCEEEEGD
## 2         143 PEPERDYLEAAIRTVIQIHMCEEEEGD PEPERDYLEAAIRAVIQIHMCEEEEGD
## 3         143 PEPERDYLEAAIRTVIQIHMCEEEEGD PEPERDYLEAAIRAVIQIHMCEEEEGD
## 4         143 PEPERDYLEAAIRTVIQIHMCEEEEGD PEPERDYLEAAIRAVIQIHMCEEEEGD
## 5         143 PEPERDYLEAAIRTVIQIHMCEEEEGD PEPERDYLEAAIRAVIQIHMCEEEEGD
## 6         143 PEPERDYLEAAIRTVIQIHMCEEEEGD PEPERDYLEAAIRAVIQIHMCEEEEGD
##   Total_RNA Tumor_RNA_Ratio Tumor_RNA Tumor_RNA_based_on_DNA MutRatio
## 1   1.35204              NA        NA              0.6576249       NA
## 2   1.35204              NA        NA              0.6576249       NA
## 3   1.35204              NA        NA              0.6576249       NA
## 4   1.35204              NA        NA              0.6576249       NA
## 5   1.35204              NA        NA              0.6576249       NA
## 6   1.35204              NA        NA              0.6576249       NA
##   MutRatio_Min MutRatio_Max
## 1           NA           NA
## 2           NA           NA
## 3           NA           NA
## 4           NA           NA
## 5           NA           NA
## 6           NA           NA
```

sample_result_INDEL_CLASS1_ALL.txt
```
##           HLA Pos      Gene Evaluated_Mutant_Peptide_Core
## 1 HLA-A*02:01   1 0_UGT2B28                    GIPMVGIPLV
## 2 HLA-A*02:01   2 0_UGT2B28                 YHGIPMVGIPLVL
## 3 HLA-A*02:01   3 0_UGT2B28                  HGIPMVGIPLVL
## 4 HLA-A*02:01   4 0_UGT2B28                    GIPMVGIPLV
## 5 HLA-A*02:01   5 0_UGT2B28                    IPMVGIPLVL
## 6 HLA-A*02:01   2 0_UGT2B28                    GIPMVGIPLV
##   Evaluated_Mutant_Peptide Mut_IC50 Mut_Rank Chr     NM_ID         Change
## 1            IYHGIPMVGIPLV   8158.7  12.6683   4 NM_053039 Out_c.1186dupT
## 2            YHGIPMVGIPLVL    264.7   1.9668   4 NM_053039 Out_c.1186dupT
## 3            HGIPMVGIPLVLG  14185.4  18.6104   4 NM_053039 Out_c.1186dupT
## 4            GIPMVGIPLVLGS  13299.6  17.6940   4 NM_053039 Out_c.1186dupT
## 5            IPMVGIPLVLGST  17674.7  22.4056   4 NM_053039 Out_c.1186dupT
## 6             YHGIPMVGIPLV   5439.5   9.9316   4 NM_053039 Out_c.1186dupT
##   Ref Alt Prob Mutation_Prob. Exon_Start Exon_End Mutation_Position
## 1   -   T    0              0   70146216 70160768        4_70156404
## 2   -   T    0              0   70146216 70160768        4_70156404
## 3   -   T    0              0   70146216 70160768        4_70156404
## 4   -   T    0              0   70146216 70160768        4_70156404
## 5   -   T    0              0   70146216 70160768        4_70156404
## 6   -   T    0              0   70146216 70160768        4_70156404
##   Total_Depth Tumor_Depth                   Wt_Peptide     Mutant_Peptide
## 1          84          43 IYHGIPMVGIPLFWDQPDNIAHMKAKGA IYHGIPMVGIPLVLGSTX
## 2          84          43 IYHGIPMVGIPLFWDQPDNIAHMKAKGA IYHGIPMVGIPLVLGSTX
## 3          84          43 IYHGIPMVGIPLFWDQPDNIAHMKAKGA IYHGIPMVGIPLVLGSTX
## 4          84          43 IYHGIPMVGIPLFWDQPDNIAHMKAKGA IYHGIPMVGIPLVLGSTX
## 5          84          43 IYHGIPMVGIPLFWDQPDNIAHMKAKGA IYHGIPMVGIPLVLGSTX
## 6          84          43 IYHGIPMVGIPLFWDQPDNIAHMKAKGA IYHGIPMVGIPLVLGSTX
##   Total_RNA Tumor_RNA_Ratio Tumor_RNA Tumor_RNA_based_on_DNA MutRatio
## 1         0              NA        NA                      0       NA
## 2         0              NA        NA                      0       NA
## 3         0              NA        NA                      0       NA
## 4         0              NA        NA                      0       NA
## 5         0              NA        NA                      0       NA
## 6         0              NA        NA                      0       NA
##   MutRatio_Min MutRatio_Max
## 1           NA           NA
## 2           NA           NA
## 3           NA           NA
## 4           NA           NA
## 5           NA           NA
## 6           NA           NA
```

sample_result_SVFusion_CLASS1_ALL.txt
```
##           HLA Pos       Gene Evaluated_Mutant_Peptide_Core
## 1 HLA-A*02:01   1 0_SLC25A12                 GDPHELRNIFLQL
## 2 HLA-A*02:01   2 0_SLC25A12                  DPHELRNIFLQL
## 3 HLA-A*02:01   3 0_SLC25A12                   ELRNIFLQLSA
## 4 HLA-A*02:01   4 0_SLC25A12                 HELRNIFLQLSAV
## 5 HLA-A*02:01   5 0_SLC25A12                  ELRNIFLQLSAV
## 6 HLA-A*02:01   6 0_SLC25A12                   LRNIFLQLSAV
##   Evaluated_Mutant_Peptide Mut_IC50 Mut_Rank Chr               NM_ID
## 1            GDPHELRNIFLQL   7994.3  12.5026   2 NM_003705_NM_015090
## 2            DPHELRNIFLQLS  36655.0  59.0000   2 NM_003705_NM_015090
## 3            PHELRNIFLQLSA  31164.0  43.8594   2 NM_003705_NM_015090
## 4            HELRNIFLQLSAV   2243.7   6.1823   2 NM_003705_NM_015090
## 5            ELRNIFLQLSAVQ  20527.4  25.8717   2 NM_003705_NM_015090
## 6            LRNIFLQLSAVQE  22787.8  28.8756   2 NM_003705_NM_015090
##                            Change Ref            Alt Prob Mutation_Prob.
## 1 In_SLC25A12_intron_NFASC_intron   A [1:204908711[A    0              0
## 2 In_SLC25A12_intron_NFASC_intron   A [1:204908711[A    0              0
## 3 In_SLC25A12_intron_NFASC_intron   A [1:204908711[A    0              0
## 4 In_SLC25A12_intron_NFASC_intron   A [1:204908711[A    0              0
## 5 In_SLC25A12_intron_NFASC_intron   A [1:204908711[A    0              0
## 6 In_SLC25A12_intron_NFASC_intron   A [1:204908711[A    0              0
##   Exon_Start  Exon_End Mutation_Position Total_Depth Tumor_Depth
## 1  172639914 172750816       2_172743385           0           0
## 2  172639914 172750816       2_172743385           0           0
## 3  172639914 172750816       2_172743385           0           0
## 4  172639914 172750816       2_172743385           0           0
## 5  172639914 172750816       2_172743385           0           0
## 6  172639914 172750816       2_172743385           0           0
##   Wt_Peptide
## 1 MAVKVQTT...
## 2 MAVKVQTT...
## 3 MAVKVQTT...
## 4 MAVKVQTT...
## 5 MAVKVQTT...
## 6 MAVKVQTT...
##                                                      Mutant_Peptide
## 1 GDPHELRNIFLQLSAVQEAQLKRLEVTRPRVLGSREQGQVPRMARQPPPPWVHAAFLLCLLSLGG
## 2 GDPHELRNIFLQLSAVQEAQLKRLEVTRPRVLGSREQGQVPRMARQPPPPWVHAAFLLCLLSLGG
## 3 GDPHELRNIFLQLSAVQEAQLKRLEVTRPRVLGSREQGQVPRMARQPPPPWVHAAFLLCLLSLGG
## 4 GDPHELRNIFLQLSAVQEAQLKRLEVTRPRVLGSREQGQVPRMARQPPPPWVHAAFLLCLLSLGG
## 5 GDPHELRNIFLQLSAVQEAQLKRLEVTRPRVLGSREQGQVPRMARQPPPPWVHAAFLLCLLSLGG
## 6 GDPHELRNIFLQLSAVQEAQLKRLEVTRPRVLGSREQGQVPRMARQPPPPWVHAAFLLCLLSLGG
##   Total_RNA Tumor_RNA_Ratio Tumor_RNA Tumor_RNA_based_on_DNA MutRatio
## 1   5.24517              NA        NA                    NaN       NA
## 2   5.24517              NA        NA                    NaN       NA
## 3   5.24517              NA        NA                    NaN       NA
## 4   5.24517              NA        NA                    NaN       NA
## 5   5.24517              NA        NA                    NaN       NA
## 6   5.24517              NA        NA                    NaN       NA
##   MutRatio_Min MutRatio_Max
## 1           NA           NA
## 2           NA           NA
## 3           NA           NA
## 4           NA           NA
## 5           NA           NA
## 6           NA           NA
```

sample_result_SeqFragment_CLASS1_ALL.txt
```
##           HLA Pos    Gene Evaluated_Mutant_Peptide_Core
## 1 HLA-A*02:01   1 0_NFKB1                 MAEDDPYLGRPEK
## 2 HLA-A*02:01   2 0_NFKB1                 AEDDPYLGRPEKM
## 3 HLA-A*02:01   3 0_NFKB1                     YLGRPEKMF
## 4 HLA-A*02:01   4 0_NFKB1                    YLGRPEKMFH
## 5 HLA-A*02:01   5 0_NFKB1                   YLGRPEKMFHL
## 6 HLA-A*02:01   6 0_NFKB1                   YLGRPEKMFHL
##   Evaluated_Mutant_Peptide Mut_IC50 Mut_Rank       Chr
## 1            MAEDDPYLGRPEK  34529.5  52.6059 chr4;chr4
## 2            AEDDPYLGRPEKM  37395.3  61.4560 chr4;chr4
## 3            EDDPYLGRPEKMF  41328.0  76.2274 chr4;chr4
## 4            DDPYLGRPEKMFH  42864.6  82.5363 chr4;chr4
## 5            DPYLGRPEKMFHL   3345.9   7.5921 chr4;chr4
## 6            PYLGRPEKMFHLD  19715.7  24.8190 chr4;chr4
##                    NM_ID      Change Ref Alt Prob Mutation_Prob.
## 1 NM_003998;NM_001165412 NFKB1;NFKB1  NA  NA   NA             NA
## 2 NM_003998;NM_001165412 NFKB1;NFKB1  NA  NA   NA             NA
## 3 NM_003998;NM_001165412 NFKB1;NFKB1  NA  NA   NA             NA
## 4 NM_003998;NM_001165412 NFKB1;NFKB1  NA  NA   NA             NA
## 5 NM_003998;NM_001165412 NFKB1;NFKB1  NA  NA   NA             NA
## 6 NM_003998;NM_001165412 NFKB1;NFKB1  NA  NA   NA             NA
##            Exon_Start            Exon_End Mutation_Position Total_Depth
## 1 103422485;103422485 103538459;103538459      chr4;chr4_NA          NA
## 2 103422485;103422485 103538459;103538459      chr4;chr4_NA          NA
## 3 103422485;103422485 103538459;103538459      chr4;chr4_NA          NA
## 4 103422485;103422485 103538459;103538459      chr4;chr4_NA          NA
## 5 103422485;103422485 103538459;103538459      chr4;chr4_NA          NA
## 6 103422485;103422485 103538459;103538459      chr4;chr4_NA          NA
##   Tumor_Depth
## 1          NA
## 2          NA
## 3          NA
## 4          NA
## 5          NA
## 6          NA
##   Wt_Peptide
## 1 MAEDDPY...
## 2 MAEDDPY...
## 3 MAEDDPY...
## 4 MAEDDPY...
## 5 MAEDDPY...
## 6 MAEDDPY...
##              Mutant_Peptide Total_RNA Tumor_RNA_Ratio Tumor_RNA
## 1 MAEDDPYLGRPEKMFHLDPSLTHTI        NA              NA        NA
## 2 MAEDDPYLGRPEKMFHLDPSLTHTI        NA              NA        NA
## 3 MAEDDPYLGRPEKMFHLDPSLTHTI        NA              NA        NA
## 4 MAEDDPYLGRPEKMFHLDPSLTHTI        NA              NA        NA
## 5 MAEDDPYLGRPEKMFHLDPSLTHTI        NA              NA        NA
## 6 MAEDDPYLGRPEKMFHLDPSLTHTI        NA              NA        NA
##   Tumor_RNA_based_on_DNA MutRatio MutRatio_Min MutRatio_Max
## 1                     NA       NA           NA           NA
## 2                     NA       NA           NA           NA
## 3                     NA       NA           NA           NA
## 4                     NA       NA           NA           NA
## 5                     NA       NA           NA           NA
## 6                     NA       NA           NA           NA
```
