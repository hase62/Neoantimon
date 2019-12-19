## A manuscript is at bioarxiv (https://www.biorxiv.org/content/10.1101/869388v1). 

## 1. Preparation
### -Install wget  (Required)
```
(CentOS) yum install wget
(Mac) brew install wget
```

### -Download and Set netMHCpan4.0 (Required)

1. Download netMHCpan4.0 from http://www.cbs.dtu.dk/cgi-bin/nph-sw_request?netMHCpan and move it to the working directory. 

2. Do initial setting at the working directory as followings.
```
wget --no-check-certificate https://github.com/hase62/Neoantimon/raw/master/lib/setNetMHCpan4.0.sh
chmod 750 setNetMHCpan4.0.sh
./setNetMHCpan4.0.sh 
```
### -Download and Set mhcflurry (Not Required)
1. (Recommended) Install anaconda from https://www.anaconda.com/distribution/, and then run the following codes. 
```
pip install mhcflurry
mhcflurry-downloads fetch
pip install mhctools
```

2. Otherwise, install python from https://www.python.org/downloads/release/python-380/, and then run the above codes.

### -Download and Set netMHCIIpan3.2 (Required)

1. Download netMHCIIpan 3.2 from http://www.cbs.dtu.dk/cgi-bin/nph-sw_request?netMHCIIpan and move it to the working directory. 

2. Do initial setting at the working directory as followings.
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
1. (Recommended) Install anaconda from https://www.anaconda.com/distribution/, and then run the following codes. 
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
#<Suggest>
install.packages('data.table');
library(data.table);
```

## 3. Data Format
### -HLA Table (or, You can indicate HLA types directly)

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

**An annovar format or an Ensembl Variant Effect Predictor (VEP) format annotated VCF file is required for Snv/Indel.**

When indicating annovar format, it must include columns representing "Chromosome Number", "Mutation Start Position", "Mutation End Position", "Mutation Ref", "Mutation Alt", and "NM_ID (AAChange.refGene)".
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
  Result_HLA1_SNV <- MainSNVClass1(input_annovar_format_file = "data/sample_vcf.txt",
                                   file_name_in_hla_table = "sample",
                                   hla_file = "data/sample_hla_table_c1.txt",
                                   refflat_file  = "lib/refFlat.txt",
                                   refmrna_file = "lib/refMrna.fa",
                                   rnaexp_file = "data/sample_rna_exp.txt",
                                   netMHCpan_dir = "lib/netMHCpan-4.0/netMHCpan",
                                   depth_tumor_column = 12,
                                   depth_normal_column = 14,
                                   SNPs = "lib/sample.snps.vcf",
                                   multiple_variants = TRUE,
                                   MHCflurry = "mhctools")

  Result_HLA2_SNV <- MainSNVClass2(input_annovar_format_file = "data/sample_vcf.txt",
                                   file_name_in_hla_table = "sample",
                                   hla_file = "data/sample_hla_table_c2.txt",
                                   refflat_file  = "lib/refFlat.txt",
                                   refmrna_file = "lib/refMrna.fa",
                                   rnaexp_file = "data/sample_rna_exp.txt",
                                   netMHCIIpan_dir = "lib/netMHCIIpan-3.2/netMHCIIpan",
                                   depth_tumor_column = 12,
                                   depth_normal_column = 14,
                                   SNPs = "lib/sample.snps.vcf",
                                   multiple_variants = TRUE)

  Result_HLA1_INDEL <- MainINDELClass1(input_annovar_format_file = "data/sample_vcf.txt",
                                       file_name_in_hla_table = "sample",
                                       hla_file = "data/sample_hla_table_c1.txt",
                                       refflat_file  = "lib/refFlat.txt",
                                       refmrna_file = "lib/refMrna.fa",
                                       rnaexp_file = "data/sample_rna_exp.txt",
                                       netMHCpan_dir = "lib/netMHCpan-4.0/netMHCpan",
                                       depth_tumor_column = 12,
                                       depth_normal_column = 14,
                                       SNPs = "lib/sample.snps.vcf",
                                       multiple_variants = TRUE,
                                       MHCflurry = "mhctools")

  Result_HLA2_INDEL <- MainINDELClass2(input_annovar_format_file = "data/sample_vcf.txt",
                                       file_name_in_hla_table = "sample",
                                       hla_file = "data/sample_hla_table_c2.txt",
                                       refflat_file  = "lib/refFlat.txt",
                                       refmrna_file = "lib/refMrna.fa",
                                       rnaexp_file = "data/sample_rna_exp.txt",
                                       netMHCIIpan_dir = "lib/netMHCIIpan-3.2/netMHCIIpan",
                                       depth_tumor_column = 12,
                                       depth_normal_column = 14,
                                       SNPs = "lib/sample.snps.vcf",
                                       multiple_variants = TRUE)
```

##### Calculate Neoantigens on SV fusions for HLA Class I and II. 
```
  Result_HLA1_SV <- MainSVFUSIONClass1(input_file = "data/sample_sv_bnd.txt",
                                       file_name_in_hla_table = "sample",
                                       hla_file = "data/sample_hla_table_c1.txt",
                                       refflat_file  = "lib/refFlat.txt",
                                       refmrna_file = "lib/refMrna.fa",
                                       rnaexp_file = "data/sample_rna_exp.txt",
                                       netMHCpan_dir = "lib/netMHCpan-4.0/netMHCpan",
                                       refdna_file = "lib/GRCh37.fa",
                                       mutation_alt_bnd_column = 5,
                                       gene_symbol_column = 7,
                                       mate_id_column = 8)

  Result_HLA2_SV <- MainSVFUSIONClass2(input_file = "data/sample_sv_bnd.txt",
                                       file_name_in_hla_table = "sample",
                                       hla_file = "data/sample_hla_table_c2.txt",
                                       refflat_file  = "lib/refFlat.txt",
                                       refmrna_file = "lib/refMrna.fa",
                                       rnaexp_file = "data/sample_rna_exp.txt",
                                       netMHCIIpan_dir = "lib/netMHCIIpan-3.2/netMHCIIpan",
                                       refdna_file = "lib/GRCh37.fa",
                                       mutation_alt_bnd_column = 5,
                                       gene_symbol_column = 7,
                                       mate_id_column = 8)
```

##### Calculate Neoantigens from a fragment of RNA sequence for HLA Class I and II by comparing to the original protein. 
```
  Result_HLA1_Seq <- MainSeqFragmentClass1(input_sequence = "atggcagaagatgatccatatttgggaaggcctgaaaaaatgtttcatttggatccttctttgactcatacaatatttaatc",
                                           file_name_in_hla_table = "sample",
                                           hla_file = "data/sample_hla_table_c1.txt",
                                           hmdir = getwd(),
                                           job_id = "NO_job_id",
                                           refflat_file  = "lib/refFlat.txt",
                                           refmrna_file = "lib/refMrna.fa",
                                           netMHCpan_dir = "lib/netMHCpan-4.0/netMHCpan",
                                           reference_nm_id = c("NM_003998", "NM_001165412"))

  Result_HLA2_Seq <- MainSeqFragmentClass2(input_sequence = "atggcagaagatgatccatatttgggaaggcctgaacaaatgtttcatttgatccttctttgactcatacaatatttaatc",
                                           file_name_in_hla_table = "sample",
                                           hla_file = "data/sample_hla_table_c2.txt",
                                           hmdir = getwd(),
                                           job_id = "NO_job_id",
                                           refflat_file  = "lib/refFlat.txt",
                                           refmrna_file = "lib/refMrna.fa",
                                           netMHCIIpan_dir = "lib/netMHCIIpan-3.2/netMHCIIpan",
                                           reference_gene_symbol = c("NFKB1", "BCL3"))
```

## 5. Result
```
print(head(Result_HLA1_SNV[[1]]))
```
  
```
##      HLA           Pos Gene        Evaluated_Mutant_Peptide Mut_IC50  Mut_Rank 
## [1,] "HLA-A*02:01" "2" "0_CYP4A11" "HQERCREEIHSLP"          "38277.9" "64.5038"
## [2,] "HLA-A*02:01" "3" "0_CYP4A11" "QERCREEIHSLPG"          "44731.5" "90.0840"
## [3,] "HLA-A*02:01" "4" "0_CYP4A11" "ERCREEIHSLPGD"          "45387.7" "92.3912"
## [4,] "HLA-A*02:01" "5" "0_CYP4A11" "RCREEIHSLPGDG"          "45916.2" "94.2252"
## [5,] "HLA-A*02:01" "6" "0_CYP4A11" "CREEIHSLPGDGA"          "41190.5" "75.6598"
## [6,] "HLA-A*02:01" "7" "0_CYP4A11" "REEIHSLPGDGAS"          "42224.3" "79.8820"
##      Evaluated_Wt_Peptide Wt_IC50   Wt_Rank   Chr NM_ID       Change     Ref
## [1,] "HQERCREEIHSLL"      "21979.5" "27.7445" "1" "NM_000778" "c.T1064C" "A"
## [2,] "QERCREEIHSLLG"      "43057.5" "83.3295" "1" "NM_000778" "c.T1064C" "A"
## [3,] "ERCREEIHSLLGD"      "43941.4" "86.9293" "1" "NM_000778" "c.T1064C" "A"
## [4,] "RCREEIHSLLGDG"      "43591.4" "85.5092" "1" "NM_000778" "c.T1064C" "A"
## [5,] "CREEIHSLLGDGA"      "40373.3" "72.4303" "1" "NM_000778" "c.T1064C" "A"
## [6,] "REEIHSLLGDGAS"      "41200.3" "75.7003" "1" "NM_000778" "c.T1064C" "A"
##      Alt Prob Mutation_Prob. Exon_Start Exon_End   Mutation_Position
## [1,] "G" "0"  "0"            "47394859" "47407148" "1_47399872"     
## [2,] "G" "0"  "0"            "47394859" "47407148" "1_47399872"     
## [3,] "G" "0"  "0"            "47394859" "47407148" "1_47399872"     
## [4,] "G" "0"  "0"            "47394859" "47407148" "1_47399872"     
## [5,] "G" "0"  "0"            "47394859" "47407148" "1_47399872"     
## [6,] "G" "0"  "0"            "47394859" "47407148" "1_47399872"     
##      Total_Depth Tumor_Depth Wt_Peptide                   
## [1,] "113"       "64"        "KHQERCREEIHSLLGDGASITWNHLDQ"
## [2,] "113"       "64"        "KHQERCREEIHSLLGDGASITWNHLDQ"
## [3,] "113"       "64"        "KHQERCREEIHSLLGDGASITWNHLDQ"
## [4,] "113"       "64"        "KHQERCREEIHSLLGDGASITWNHLDQ"
## [5,] "113"       "64"        "KHQERCREEIHSLLGDGASITWNHLDQ"
## [6,] "113"       "64"        "KHQERCREEIHSLLGDGASITWNHLDQ"
##      Mutant_Peptide                Total_RNA Tumor_RNA_Ratio Tumor_RNA
## [1,] "KHQERCREEIHSLPGDGASITWNHLDQ" "0"       "NA"            "NA"     
## [2,] "KHQERCREEIHSLPGDGASITWNHLDQ" "0"       "NA"            "NA"     
## [3,] "KHQERCREEIHSLPGDGASITWNHLDQ" "0"       "NA"            "NA"     
## [4,] "KHQERCREEIHSLPGDGASITWNHLDQ" "0"       "NA"            "NA"     
## [5,] "KHQERCREEIHSLPGDGASITWNHLDQ" "0"       "NA"            "NA"     
## [6,] "KHQERCREEIHSLPGDGASITWNHLDQ" "0"       "NA"            "NA"     
##      Tumor_RNA_based_on_DNA MutRatio MutRatio_Min MutRatio_Max
## [1,] "0"                    "NA"     "NA"         "NA"        
## [2,] "0"                    "NA"     "NA"         "NA"        
## [3,] "0"                    "NA"     "NA"         "NA"        
## [4,] "0"                    "NA"     "NA"         "NA"        
## [5,] "0"                    "NA"     "NA"         "NA"        
## [6,] "0"                    "NA"     "NA"         "NA"
```

```
print(Export_Summary_SNV(Input = Result_HLA1_SNV[[1]], Mut_IC50_th = 500, Wt_IC50_th = 500))
```

```
##              Num_All_Alteration        Num_Evaluated_Alteration 
##                               4                               4 
## Num_Alteration_Generating_NeoAg                 Num_All_Peptide 
##                               2                             252 
##           Num_Evaluated_Peptide    Num_Peptide_Generating_NeoAg 
##                             252                               3
```

```
print(head(Result_HLA2_SNV))
```
  
```
##      HLA                     Pos Gene        Evaluated_Mutant_Peptide
## [1,] "HLA-DPA10103-DPB10201" "2" "0_CYP4A11" "PKHQERCREEIHSLP"       
## [2,] "HLA-DPA10103-DPB10201" "3" "0_CYP4A11" "KHQERCREEIHSLPG"       
## [3,] "HLA-DPA10103-DPB10201" "4" "0_CYP4A11" "HQERCREEIHSLPGD"       
## [4,] "HLA-DPA10103-DPB10201" "5" "0_CYP4A11" "QERCREEIHSLPGDG"       
## [5,] "HLA-DPA10103-DPB10201" "6" "0_CYP4A11" "ERCREEIHSLPGDGA"       
## [6,] "HLA-DPA10103-DPB10201" "7" "0_CYP4A11" "RCREEIHSLPGDGAS"       
##      Mut_IC50   Mut_Rank Evaluated_Wt_Peptide Wt_IC50   Wt_Rank Chr NM_ID      
## [1,] "12151.54" "95.00"  "PKHQERCREEIHSLL"    "5209.84" "75.00" "1" "NM_000778"
## [2,] "11550.05" "90.00"  "KHQERCREEIHSLLG"    "4635.67" "70.00" "1" "NM_000778"
## [3,] "12677.06" "95.00"  "HQERCREEIHSLLGD"    "4880.31" "75.00" "1" "NM_000778"
## [4,] "13675.38" "95.00"  "QERCREEIHSLLGDG"    "5435.50" "75.00" "1" "NM_000778"
## [5,] "13240.19" "95.00"  "ERCREEIHSLLGDGA"    "5173.27" "75.00" "1" "NM_000778"
## [6,] "13790.24" "95.00"  "RCREEIHSLLGDGAS"    "5387.49" "75.00" "1" "NM_000778"
##      Change     Ref Alt Prob Mutation_Prob. Exon_Start Exon_End  
## [1,] "c.T1064C" "A" "G" "0"  "0"            "47394859" "47407148"
## [2,] "c.T1064C" "A" "G" "0"  "0"            "47394859" "47407148"
## [3,] "c.T1064C" "A" "G" "0"  "0"            "47394859" "47407148"
## [4,] "c.T1064C" "A" "G" "0"  "0"            "47394859" "47407148"
## [5,] "c.T1064C" "A" "G" "0"  "0"            "47394859" "47407148"
## [6,] "c.T1064C" "A" "G" "0"  "0"            "47394859" "47407148"
##      Mutation_Position Total_Depth Tumor_Depth
## [1,] "1_47399872"      "113"       "64"       
## [2,] "1_47399872"      "113"       "64"       
## [3,] "1_47399872"      "113"       "64"       
## [4,] "1_47399872"      "113"       "64"       
## [5,] "1_47399872"      "113"       "64"       
## [6,] "1_47399872"      "113"       "64"       
##      Wt_Peptide                        Mutant_Peptide                   
## [1,] "HPKHQERCREEIHSLLGDGASITWNHLDQMP" "HPKHQERCREEIHSLPGDGASITWNHLDQMP"
## [2,] "HPKHQERCREEIHSLLGDGASITWNHLDQMP" "HPKHQERCREEIHSLPGDGASITWNHLDQMP"
## [3,] "HPKHQERCREEIHSLLGDGASITWNHLDQMP" "HPKHQERCREEIHSLPGDGASITWNHLDQMP"
## [4,] "HPKHQERCREEIHSLLGDGASITWNHLDQMP" "HPKHQERCREEIHSLPGDGASITWNHLDQMP"
## [5,] "HPKHQERCREEIHSLLGDGASITWNHLDQMP" "HPKHQERCREEIHSLPGDGASITWNHLDQMP"
## [6,] "HPKHQERCREEIHSLLGDGASITWNHLDQMP" "HPKHQERCREEIHSLPGDGASITWNHLDQMP"
##      Total_RNA Tumor_RNA_Ratio Tumor_RNA Tumor_RNA_based_on_DNA MutRatio
## [1,] "0"       "NA"            "NA"      "0"                    "NA"    
## [2,] "0"       "NA"            "NA"      "0"                    "NA"    
## [3,] "0"       "NA"            "NA"      "0"                    "NA"    
## [4,] "0"       "NA"            "NA"      "0"                    "NA"    
## [5,] "0"       "NA"            "NA"      "0"                    "NA"    
## [6,] "0"       "NA"            "NA"      "0"                    "NA"    
##      MutRatio_Min MutRatio_Max
## [1,] "NA"         "NA"        
## [2,] "NA"         "NA"        
## [3,] "NA"         "NA"        
## [4,] "NA"         "NA"        
## [5,] "NA"         "NA"        
## [6,] "NA"         "NA"
```

```
print(Export_Summary_SNV(Input = Result_HLA2_SNV, Mut_IC50_th = 500, Wt_IC50_th = 500))
```

```
##              Num_All_Alteration        Num_Evaluated_Alteration 
##                               4                               4 
## Num_Alteration_Generating_NeoAg                 Num_All_Peptide 
##                               2                              60 
##           Num_Evaluated_Peptide    Num_Peptide_Generating_NeoAg 
##                              60                              10
```

```
print(head(Result_HLA1_INDEL[[1]]))
```

```
##      HLA           Pos Gene        Evaluated_Mutant_Peptide_Core
## [1,] "HLA-A*02:01" "1" "0_UGT2B28" "GIPMVGIPLV"                 
## [2,] "HLA-A*02:01" "2" "0_UGT2B28" "YHGIPMVGIPLVL"              
## [3,] "HLA-A*02:01" "3" "0_UGT2B28" "HGIPMVGIPLVL"               
## [4,] "HLA-A*02:01" "4" "0_UGT2B28" "GIPMVGIPLV"                 
## [5,] "HLA-A*02:01" "5" "0_UGT2B28" "IPMVGIPLVL"                 
## [6,] "HLA-A*02:01" "2" "0_UGT2B28" "GIPMVGIPLV"                 
##      Evaluated_Mutant_Peptide Mut_IC50  Mut_Rank  Chr NM_ID      
## [1,] "IYHGIPMVGIPLV"          "8158.7"  "12.6683" "4" "NM_053039"
## [2,] "YHGIPMVGIPLVL"          "264.7"   "1.9668"  "4" "NM_053039"
## [3,] "HGIPMVGIPLVLG"          "14185.4" "18.6104" "4" "NM_053039"
## [4,] "GIPMVGIPLVLGS"          "13299.6" "17.6940" "4" "NM_053039"
## [5,] "IPMVGIPLVLGST"          "17674.7" "22.4056" "4" "NM_053039"
## [6,] "YHGIPMVGIPLV"           "5439.5"  "9.9316"  "4" "NM_053039"
##      Change           Ref Alt Prob Mutation_Prob. Exon_Start Exon_End  
## [1,] "Out_c.1186dupT" "-" "T" "0"  "0"            "70146192" "70160768"
## [2,] "Out_c.1186dupT" "-" "T" "0"  "0"            "70146192" "70160768"
## [3,] "Out_c.1186dupT" "-" "T" "0"  "0"            "70146192" "70160768"
## [4,] "Out_c.1186dupT" "-" "T" "0"  "0"            "70146192" "70160768"
## [5,] "Out_c.1186dupT" "-" "T" "0"  "0"            "70146192" "70160768"
## [6,] "Out_c.1186dupT" "-" "T" "0"  "0"            "70146192" "70160768"
##      Mutation_Position Total_Depth Tumor_Depth Wt_Peptide                    
## [1,] "4_70156404"      "84"        "43"        "IYHGIPMVGIPLFWDQPDNIAHMKAKGA"
## [2,] "4_70156404"      "84"        "43"        "IYHGIPMVGIPLFWDQPDNIAHMKAKGA"
## [3,] "4_70156404"      "84"        "43"        "IYHGIPMVGIPLFWDQPDNIAHMKAKGA"
## [4,] "4_70156404"      "84"        "43"        "IYHGIPMVGIPLFWDQPDNIAHMKAKGA"
## [5,] "4_70156404"      "84"        "43"        "IYHGIPMVGIPLFWDQPDNIAHMKAKGA"
## [6,] "4_70156404"      "84"        "43"        "IYHGIPMVGIPLFWDQPDNIAHMKAKGA"
##      Mutant_Peptide       Total_RNA Tumor_RNA_Ratio Tumor_RNA
## [1,] "IYHGIPMVGIPLVLGSTX" "0"       "NA"            "NA"     
## [2,] "IYHGIPMVGIPLVLGSTX" "0"       "NA"            "NA"     
## [3,] "IYHGIPMVGIPLVLGSTX" "0"       "NA"            "NA"     
## [4,] "IYHGIPMVGIPLVLGSTX" "0"       "NA"            "NA"     
## [5,] "IYHGIPMVGIPLVLGSTX" "0"       "NA"            "NA"     
## [6,] "IYHGIPMVGIPLVLGSTX" "0"       "NA"            "NA"     
##      Tumor_RNA_based_on_DNA MutRatio MutRatio_Min MutRatio_Max
## [1,] "0"                    "NA"     "NA"         "NA"        
## [2,] "0"                    "NA"     "NA"         "NA"        
## [3,] "0"                    "NA"     "NA"         "NA"        
## [4,] "0"                    "NA"     "NA"         "NA"        
## [5,] "0"                    "NA"     "NA"         "NA"        
## [6,] "0"                    "NA"     "NA"         "NA"
```

```
print(Export_Summary_IndelSV(Input = Result_HLA1_INDEL[[1]], Mut_IC50_th = 500))
```

```
##              Num_All_Alteration        Num_Evaluated_Alteration 
##                               3                               3 
## Num_Alteration_Generating_NeoAg                 Num_All_Peptide 
##                               3                             270 
##           Num_Evaluated_Peptide    Num_Peptide_Generating_NeoAg 
##                             270                              15
```

```
print(Export_Summary_IndelSV_perFragments(Input = Result_HLA1_INDEL[[1]], Mut_IC50_th = 500))
```

```
##                                                     IYHGIPMVGIPLVLGSTX-0_UGT2B28
## Num_Peptide_Per_Pep                                                       30.000
## Num_Cond_Peptide_Per_Pep                                                  30.000
## Num_Rest_Peptide_Per_Pep                                                   4.000
## Num_Rest_Peptide_Per_Pep / Num_Cond_Peptide_Per_Pep                        0.133
## -logP                                                                      0.414
##                                                     YQRDELLAAIKKFHIKVATQX-1_COL12A1
## Num_Peptide_Per_Pep                                                          48.000
## Num_Cond_Peptide_Per_Pep                                                     48.000
## Num_Rest_Peptide_Per_Pep                                                      2.000
## Num_Rest_Peptide_Per_Pep / Num_Cond_Peptide_Per_Pep                           0.042
## -logP                                                                         0.483
##                                                     EMARDFLPSLKNPFWKPSILPIFMYKHCSVQFSVRHGDVQTKVHX-2_SLCO1C1
## Num_Peptide_Per_Pep                                                                                 192.000
## Num_Cond_Peptide_Per_Pep                                                                            192.000
## Num_Rest_Peptide_Per_Pep                                                                              9.000
## Num_Rest_Peptide_Per_Pep / Num_Cond_Peptide_Per_Pep                                                   0.047
## -logP                                                                                                 1.035
```

```
print(head(Result_HLA2_INDEL))
```

```
##      HLA                     Pos Gene        Evaluated_Mutant_Peptide_Core
## [1,] "HLA-DPA10103-DPB10201" "3" "0_UGT2B28" "YHGIPMVGI"                  
## [2,] "HLA-DPA10103-DPB10201" "2" "0_UGT2B28" "YHGIPMVGI"                  
## [3,] "HLA-DPA10103-DPB10201" "1" "0_UGT2B28" "YHGIPMVGI"                  
## [4,] "HLA-DPA10103-DPB10201" "4" "0_UGT2B28" "PMVGIPLVL"                  
## [5,] "HLA-DPA10103-DPB10201" "3" "0_UGT2B28" "PMVGIPLVL"                  
## [6,] "HLA-DPA10103-DPB10201" "2" "1_COL12A1" "YQRDELLAA"                  
##      Evaluated_Mutant_Peptide Mut_IC50  Mut_Rank Chr NM_ID      
## [1,] "EAIYHGIPMVGIPLV"        "682.35"  "28.00"  "4" "NM_053039"
## [2,] "AIYHGIPMVGIPLVL"        "643.18"  "27.00"  "4" "NM_053039"
## [3,] "IYHGIPMVGIPLVLG"        "799.71"  "31.00"  "4" "NM_053039"
## [4,] "YHGIPMVGIPLVLGS"        "984.87"  "35.00"  "4" "NM_053039"
## [5,] "HGIPMVGIPLVLGST"        "1262.77" "40.00"  "4" "NM_053039"
## [6,] "QYYQRDELLAAIKKF"        "528.36"  "23.00"  "6" "NM_004370"
##      Change           Ref Alt Prob Mutation_Prob. Exon_Start Exon_End  
## [1,] "Out_c.1186dupT" "-" "T" "0"  "0"            "70146192" "70160768"
## [2,] "Out_c.1186dupT" "-" "T" "0"  "0"            "70146192" "70160768"
## [3,] "Out_c.1186dupT" "-" "T" "0"  "0"            "70146192" "70160768"
## [4,] "Out_c.1186dupT" "-" "T" "0"  "0"            "70146192" "70160768"
## [5,] "Out_c.1186dupT" "-" "T" "0"  "0"            "70146192" "70160768"
## [6,] "Out_c.628delA"  "T" "-" "0"  "0"            "75794041" "75915769"
##      Mutation_Position Total_Depth Tumor_Depth
## [1,] "4_70156404"      "84"        "43"       
## [2,] "4_70156404"      "84"        "43"       
## [3,] "4_70156404"      "84"        "43"       
## [4,] "4_70156404"      "84"        "43"       
## [5,] "4_70156404"      "84"        "43"       
## [6,] "6_75899298"      "195"       "122"      
##      Wt_Peptide                          Mutant_Peptide            Total_RNA
## [1,] "EAIYHGIPMVGIPLFWDQPDNIAHMKAKGA"    "EAIYHGIPMVGIPLVLGSTX"    "0"      
## [2,] "EAIYHGIPMVGIPLFWDQPDNIAHMKAKGA"    "EAIYHGIPMVGIPLVLGSTX"    "0"      
## [3,] "EAIYHGIPMVGIPLFWDQPDNIAHMKAKGA"    "EAIYHGIPMVGIPLVLGSTX"    "0"      
## [4,] "EAIYHGIPMVGIPLFWDQPDNIAHMKAKGA"    "EAIYHGIPMVGIPLVLGSTX"    "0"      
## [5,] "EAIYHGIPMVGIPLFWDQPDNIAHMKAKGA"    "EAIYHGIPMVGIPLVLGSTX"    "0"      
## [6,] "QYYQRDELLAAIKKIPYKGGNTMTGDAIDYLVK" "QYYQRDELLAAIKKFHIKVATQX" "0"      
##      Tumor_RNA_Ratio Tumor_RNA Tumor_RNA_based_on_DNA MutRatio MutRatio_Min
## [1,] "NA"            "NA"      "0"                    "NA"     "NA"        
## [2,] "NA"            "NA"      "0"                    "NA"     "NA"        
## [3,] "NA"            "NA"      "0"                    "NA"     "NA"        
## [4,] "NA"            "NA"      "0"                    "NA"     "NA"        
## [5,] "NA"            "NA"      "0"                    "NA"     "NA"        
## [6,] "NA"            "NA"      "0"                    "NA"     "NA"        
##      MutRatio_Max
## [1,] "NA"        
## [2,] "NA"        
## [3,] "NA"        
## [4,] "NA"        
## [5,] "NA"        
## [6,] "NA"
```

```
print(Export_Summary_IndelSV(Input = Result_HLA2_INDEL, Mut_IC50_th = 500))
```

```
##              Num_All_Alteration        Num_Evaluated_Alteration 
##                               3                               3 
## Num_Alteration_Generating_NeoAg                 Num_All_Peptide 
##                               3                              45 
##           Num_Evaluated_Peptide    Num_Peptide_Generating_NeoAg 
##                              45                              32
```

```
print(Export_Summary_IndelSV_perFragments(Input = Result_HLA2_INDEL, Mut_IC50_th = 500))
```

```
##                                                     EAIYHGIPMVGIPLVLGSTX-0_UGT2B28
## Num_Peptide_Per_Pep                                                           5.00
## Num_Cond_Peptide_Per_Pep                                                      5.00
## Num_Rest_Peptide_Per_Pep                                                      5.00
## Num_Rest_Peptide_Per_Pep / Num_Cond_Peptide_Per_Pep                           1.00
## -logP                                                                         0.46
##                                                     QYYQRDELLAAIKKFHIKVATQX-1_COL12A1
## Num_Peptide_Per_Pep                                                             8.000
## Num_Cond_Peptide_Per_Pep                                                        8.000
## Num_Rest_Peptide_Per_Pep                                                        7.000
## Num_Rest_Peptide_Per_Pep / Num_Cond_Peptide_Per_Pep                             0.875
## -logP                                                                           0.529
##                                                     IMEMARDFLPSLKNPFWKPSILPIFMYKHCSVQFSVRHGDVQTKVHX-2_SLCO1C1
## Num_Peptide_Per_Pep                                                                                    32.000
## Num_Cond_Peptide_Per_Pep                                                                               32.000
## Num_Rest_Peptide_Per_Pep                                                                               20.000
## Num_Rest_Peptide_Per_Pep / Num_Cond_Peptide_Per_Pep                                                     0.625
## -logP                                                                                                   1.081
```

```
print(head(Result_HLA1_SV))
```
  
```
##      HLA           Pos Gene     Evaluated_Mutant_Peptide_Core
## [1,] "HLA-A*02:01" "1" "0_AAR2" "YNSWEVGPKF"                 
## [2,] "HLA-A*02:01" "2" "0_AAR2" "YNSWEVGPKF"                 
## [3,] "HLA-A*02:01" "3" "0_AAR2" "NSWEVGPKFREQL"              
## [4,] "HLA-A*02:01" "4" "0_AAR2" "WEVGPKFREQLK"               
## [5,] "HLA-A*02:01" "5" "0_AAR2" "WEVGPKFREQLKL"              
## [6,] "HLA-A*02:01" "6" "0_AAR2" "EVGPKFREQLKL"               
##      Evaluated_Mutant_Peptide Mut_IC50  Mut_Rank  Chr  NM_ID                
## [1,] "DYNSWEVGPKFRE"          "41548.2" "77.1323" "20" "NM_015511_NM_015906"
## [2,] "YNSWEVGPKFREQ"          "37159.3" "60.6291" "20" "NM_015511_NM_015906"
## [3,] "NSWEVGPKFREQL"          "17541.5" "22.2465" "20" "NM_015511_NM_015906"
## [4,] "SWEVGPKFREQLK"          "41225.3" "75.8035" "20" "NM_015511_NM_015906"
## [5,] "WEVGPKFREQLKL"          "30625.8" "42.6216" "20" "NM_015511_NM_015906"
## [6,] "EVGPKFREQLKLF"          "40083.0" "71.2831" "20" "NM_015511_NM_015906"
##      Change                     Ref Alt              Prob Mutation_Prob.
## [1,] "In_AAR2_exon_TRIM33_exon" "G" "G]1:115005805]" "0"  "0"           
## [2,] "In_AAR2_exon_TRIM33_exon" "G" "G]1:115005805]" "0"  "0"           
## [3,] "In_AAR2_exon_TRIM33_exon" "G" "G]1:115005805]" "0"  "0"           
## [4,] "In_AAR2_exon_TRIM33_exon" "G" "G]1:115005805]" "0"  "0"           
## [5,] "In_AAR2_exon_TRIM33_exon" "G" "G]1:115005805]" "0"  "0"           
## [6,] "In_AAR2_exon_TRIM33_exon" "G" "G]1:115005805]" "0"  "0"           
##      Exon_Start Exon_End   Mutation_Position Total_Depth Tumor_Depth
## [1,] "34824338" "34844863" "20_34827929"     "0"         "0"        
## [2,] "34824338" "34844863" "20_34827929"     "0"         "0"        
## [3,] "34824338" "34844863" "20_34827929"     "0"         "0"        
## [4,] "34824338" "34844863" "20_34827929"     "0"         "0"        
## [5,] "34824338" "34844863" "20_34827929"     "0"         "0"        
## [6,] "34824338" "34844863" "20_34827929"     "0"         "0"        
##      Wt_Peptide                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  
## [1,] "MAAVQMDPELAKRLFFEGATVVILNMPKGTEFGIDYNSWEVGPKFRGVKMIPPGIHFLHYSSVDKANPKEVGPRMGFFLSLHQRGLTVLRWSTLREEVDLSPAPESEVEAMRANLQELDQFLGPYPYATLKKWISLTNFISEATVEKLQPENRQICAFSDVLPVLSMKHTKDRVGQNLPRCGIECKSYQEGLARLPEMKPRAGTEIRFSELPTQMFPEGATPAEITKHSMDLSYALETVLNKQFPSSPQDVLGELQFAFVCFLLGNVYEAFEHWKRLLNLLCRSEAAMMKHHTLYINLISILYHQLGEIPADFFVDIVSQDNFLTSTLQVFFSSACSIAVDATLRKKAEKFQAHLTKKFRWDFAAEPEDCAPVVVELPEGIEMGXXMAENKGGGEAESGGGGSGSAPVTAGAAGPAAQEAEPPLTAVLVEEEEEEGGRAGAEGGAAGPDDGGVAAASSGSAQAASSPAASVGTGVAGGAVSTPAPAPASAPAPGPSAGPPPGPPASLLDTCAVCQQSLQSRREAEPKLLPCLHSFCLRCLPEPERQLSVPIPGGSNGDIQQVGVIRCPVCRQECRQIDLVDNYFVKDTSEAPSSSDEKSEQVCTSCEDNASAVGFCVECGEWLCKTCIEAHQRVKFTKDHLIRKKEDVSESVGASGQRPVFCPVHKQEQLKLFCETCDRLTCRDCQLLEHKEHRYQFLEEAFQNQKGAIENLLAKLLEKKNYVHFAATQVQNRIKEVNETNKRVEQEIKVAIFTLINEINKKGKSLLQQLENVTKERQMKLLQQQNDITGLSRQVKHVMNFTNWAIASGSSTALLYSKRLITFQLRHILKARCDPVPAANGAIRFHCDPTFWAKNVVNLGNLVIESKPAPGYTPNVVVGQVPPGTNHISKTPGQINLAQLRLQHMQQQVYAQKHQQLQQMRMQQPPAPVPTTTTTTQQHPRQAAPQMLQQQPPRLISVQTMQRGNMNCGAFQAHQMRLAQNAARIPGIPRHSGPQYSMMQPHLQRQHSNPGHAGPFPVVSVHNTTINPTSPTTATMANANRGPTSPSVTAIELIPSVTNPENLPSLPDIPPIQLEDAGSSSLDNLLSRYISGSHLPPQPTSTMNPSPGPSALSPGSSGLSNSHTPVRPPSTSSTGSRGSCGSSGRTAEKTSLSFKSDQVKVKQEPGTEDEICSFSGGVKQEKTEDGRRSACMLSSPESSLTPPLSTNLHLESELDALASLENHVKIEPADMNESCKQSGLSSLVNGKSPIRSLMHRSARIGGDGNNKDDDPNEDWCAVCQNGGDLLCCEKCPKVFHLTCHVPTLLSFPSGDWICTFCRDIGKPEVEYDCDNLQHSKKGKTAQGLSPVDQRKCERLLLYLYCHELSIEFQEPVPASIPNYYKIIKKPMDLSTVKKKLQKKHSQHYQIPDDFVADVRLIFKNCERFNEMMKVVQVYADTQEINLKADSEVAQAGKAVALYFEDKLTEIYSDRTFAPLPEFEQEEDDGEVTEDSDEDFIQPRRKRLKSDERPVHIKX"
## [2,] "MAAVQMDPELAKRLFFEGATVVILNMPKGTEFGIDYNSWEVGPKFRGVKMIPPGIHFLHYSSVDKANPKEVGPRMGFFLSLHQRGLTVLRWSTLREEVDLSPAPESEVEAMRANLQELDQFLGPYPYATLKKWISLTNFISEATVEKLQPENRQICAFSDVLPVLSMKHTKDRVGQNLPRCGIECKSYQEGLARLPEMKPRAGTEIRFSELPTQMFPEGATPAEITKHSMDLSYALETVLNKQFPSSPQDVLGELQFAFVCFLLGNVYEAFEHWKRLLNLLCRSEAAMMKHHTLYINLISILYHQLGEIPADFFVDIVSQDNFLTSTLQVFFSSACSIAVDATLRKKAEKFQAHLTKKFRWDFAAEPEDCAPVVVELPEGIEMGXXMAENKGGGEAESGGGGSGSAPVTAGAAGPAAQEAEPPLTAVLVEEEEEEGGRAGAEGGAAGPDDGGVAAASSGSAQAASSPAASVGTGVAGGAVSTPAPAPASAPAPGPSAGPPPGPPASLLDTCAVCQQSLQSRREAEPKLLPCLHSFCLRCLPEPERQLSVPIPGGSNGDIQQVGVIRCPVCRQECRQIDLVDNYFVKDTSEAPSSSDEKSEQVCTSCEDNASAVGFCVECGEWLCKTCIEAHQRVKFTKDHLIRKKEDVSESVGASGQRPVFCPVHKQEQLKLFCETCDRLTCRDCQLLEHKEHRYQFLEEAFQNQKGAIENLLAKLLEKKNYVHFAATQVQNRIKEVNETNKRVEQEIKVAIFTLINEINKKGKSLLQQLENVTKERQMKLLQQQNDITGLSRQVKHVMNFTNWAIASGSSTALLYSKRLITFQLRHILKARCDPVPAANGAIRFHCDPTFWAKNVVNLGNLVIESKPAPGYTPNVVVGQVPPGTNHISKTPGQINLAQLRLQHMQQQVYAQKHQQLQQMRMQQPPAPVPTTTTTTQQHPRQAAPQMLQQQPPRLISVQTMQRGNMNCGAFQAHQMRLAQNAARIPGIPRHSGPQYSMMQPHLQRQHSNPGHAGPFPVVSVHNTTINPTSPTTATMANANRGPTSPSVTAIELIPSVTNPENLPSLPDIPPIQLEDAGSSSLDNLLSRYISGSHLPPQPTSTMNPSPGPSALSPGSSGLSNSHTPVRPPSTSSTGSRGSCGSSGRTAEKTSLSFKSDQVKVKQEPGTEDEICSFSGGVKQEKTEDGRRSACMLSSPESSLTPPLSTNLHLESELDALASLENHVKIEPADMNESCKQSGLSSLVNGKSPIRSLMHRSARIGGDGNNKDDDPNEDWCAVCQNGGDLLCCEKCPKVFHLTCHVPTLLSFPSGDWICTFCRDIGKPEVEYDCDNLQHSKKGKTAQGLSPVDQRKCERLLLYLYCHELSIEFQEPVPASIPNYYKIIKKPMDLSTVKKKLQKKHSQHYQIPDDFVADVRLIFKNCERFNEMMKVVQVYADTQEINLKADSEVAQAGKAVALYFEDKLTEIYSDRTFAPLPEFEQEEDDGEVTEDSDEDFIQPRRKRLKSDERPVHIKX"
## [3,] "MAAVQMDPELAKRLFFEGATVVILNMPKGTEFGIDYNSWEVGPKFRGVKMIPPGIHFLHYSSVDKANPKEVGPRMGFFLSLHQRGLTVLRWSTLREEVDLSPAPESEVEAMRANLQELDQFLGPYPYATLKKWISLTNFISEATVEKLQPENRQICAFSDVLPVLSMKHTKDRVGQNLPRCGIECKSYQEGLARLPEMKPRAGTEIRFSELPTQMFPEGATPAEITKHSMDLSYALETVLNKQFPSSPQDVLGELQFAFVCFLLGNVYEAFEHWKRLLNLLCRSEAAMMKHHTLYINLISILYHQLGEIPADFFVDIVSQDNFLTSTLQVFFSSACSIAVDATLRKKAEKFQAHLTKKFRWDFAAEPEDCAPVVVELPEGIEMGXXMAENKGGGEAESGGGGSGSAPVTAGAAGPAAQEAEPPLTAVLVEEEEEEGGRAGAEGGAAGPDDGGVAAASSGSAQAASSPAASVGTGVAGGAVSTPAPAPASAPAPGPSAGPPPGPPASLLDTCAVCQQSLQSRREAEPKLLPCLHSFCLRCLPEPERQLSVPIPGGSNGDIQQVGVIRCPVCRQECRQIDLVDNYFVKDTSEAPSSSDEKSEQVCTSCEDNASAVGFCVECGEWLCKTCIEAHQRVKFTKDHLIRKKEDVSESVGASGQRPVFCPVHKQEQLKLFCETCDRLTCRDCQLLEHKEHRYQFLEEAFQNQKGAIENLLAKLLEKKNYVHFAATQVQNRIKEVNETNKRVEQEIKVAIFTLINEINKKGKSLLQQLENVTKERQMKLLQQQNDITGLSRQVKHVMNFTNWAIASGSSTALLYSKRLITFQLRHILKARCDPVPAANGAIRFHCDPTFWAKNVVNLGNLVIESKPAPGYTPNVVVGQVPPGTNHISKTPGQINLAQLRLQHMQQQVYAQKHQQLQQMRMQQPPAPVPTTTTTTQQHPRQAAPQMLQQQPPRLISVQTMQRGNMNCGAFQAHQMRLAQNAARIPGIPRHSGPQYSMMQPHLQRQHSNPGHAGPFPVVSVHNTTINPTSPTTATMANANRGPTSPSVTAIELIPSVTNPENLPSLPDIPPIQLEDAGSSSLDNLLSRYISGSHLPPQPTSTMNPSPGPSALSPGSSGLSNSHTPVRPPSTSSTGSRGSCGSSGRTAEKTSLSFKSDQVKVKQEPGTEDEICSFSGGVKQEKTEDGRRSACMLSSPESSLTPPLSTNLHLESELDALASLENHVKIEPADMNESCKQSGLSSLVNGKSPIRSLMHRSARIGGDGNNKDDDPNEDWCAVCQNGGDLLCCEKCPKVFHLTCHVPTLLSFPSGDWICTFCRDIGKPEVEYDCDNLQHSKKGKTAQGLSPVDQRKCERLLLYLYCHELSIEFQEPVPASIPNYYKIIKKPMDLSTVKKKLQKKHSQHYQIPDDFVADVRLIFKNCERFNEMMKVVQVYADTQEINLKADSEVAQAGKAVALYFEDKLTEIYSDRTFAPLPEFEQEEDDGEVTEDSDEDFIQPRRKRLKSDERPVHIKX"
## [4,] "MAAVQMDPELAKRLFFEGATVVILNMPKGTEFGIDYNSWEVGPKFRGVKMIPPGIHFLHYSSVDKANPKEVGPRMGFFLSLHQRGLTVLRWSTLREEVDLSPAPESEVEAMRANLQELDQFLGPYPYATLKKWISLTNFISEATVEKLQPENRQICAFSDVLPVLSMKHTKDRVGQNLPRCGIECKSYQEGLARLPEMKPRAGTEIRFSELPTQMFPEGATPAEITKHSMDLSYALETVLNKQFPSSPQDVLGELQFAFVCFLLGNVYEAFEHWKRLLNLLCRSEAAMMKHHTLYINLISILYHQLGEIPADFFVDIVSQDNFLTSTLQVFFSSACSIAVDATLRKKAEKFQAHLTKKFRWDFAAEPEDCAPVVVELPEGIEMGXXMAENKGGGEAESGGGGSGSAPVTAGAAGPAAQEAEPPLTAVLVEEEEEEGGRAGAEGGAAGPDDGGVAAASSGSAQAASSPAASVGTGVAGGAVSTPAPAPASAPAPGPSAGPPPGPPASLLDTCAVCQQSLQSRREAEPKLLPCLHSFCLRCLPEPERQLSVPIPGGSNGDIQQVGVIRCPVCRQECRQIDLVDNYFVKDTSEAPSSSDEKSEQVCTSCEDNASAVGFCVECGEWLCKTCIEAHQRVKFTKDHLIRKKEDVSESVGASGQRPVFCPVHKQEQLKLFCETCDRLTCRDCQLLEHKEHRYQFLEEAFQNQKGAIENLLAKLLEKKNYVHFAATQVQNRIKEVNETNKRVEQEIKVAIFTLINEINKKGKSLLQQLENVTKERQMKLLQQQNDITGLSRQVKHVMNFTNWAIASGSSTALLYSKRLITFQLRHILKARCDPVPAANGAIRFHCDPTFWAKNVVNLGNLVIESKPAPGYTPNVVVGQVPPGTNHISKTPGQINLAQLRLQHMQQQVYAQKHQQLQQMRMQQPPAPVPTTTTTTQQHPRQAAPQMLQQQPPRLISVQTMQRGNMNCGAFQAHQMRLAQNAARIPGIPRHSGPQYSMMQPHLQRQHSNPGHAGPFPVVSVHNTTINPTSPTTATMANANRGPTSPSVTAIELIPSVTNPENLPSLPDIPPIQLEDAGSSSLDNLLSRYISGSHLPPQPTSTMNPSPGPSALSPGSSGLSNSHTPVRPPSTSSTGSRGSCGSSGRTAEKTSLSFKSDQVKVKQEPGTEDEICSFSGGVKQEKTEDGRRSACMLSSPESSLTPPLSTNLHLESELDALASLENHVKIEPADMNESCKQSGLSSLVNGKSPIRSLMHRSARIGGDGNNKDDDPNEDWCAVCQNGGDLLCCEKCPKVFHLTCHVPTLLSFPSGDWICTFCRDIGKPEVEYDCDNLQHSKKGKTAQGLSPVDQRKCERLLLYLYCHELSIEFQEPVPASIPNYYKIIKKPMDLSTVKKKLQKKHSQHYQIPDDFVADVRLIFKNCERFNEMMKVVQVYADTQEINLKADSEVAQAGKAVALYFEDKLTEIYSDRTFAPLPEFEQEEDDGEVTEDSDEDFIQPRRKRLKSDERPVHIKX"
## [5,] "MAAVQMDPELAKRLFFEGATVVILNMPKGTEFGIDYNSWEVGPKFRGVKMIPPGIHFLHYSSVDKANPKEVGPRMGFFLSLHQRGLTVLRWSTLREEVDLSPAPESEVEAMRANLQELDQFLGPYPYATLKKWISLTNFISEATVEKLQPENRQICAFSDVLPVLSMKHTKDRVGQNLPRCGIECKSYQEGLARLPEMKPRAGTEIRFSELPTQMFPEGATPAEITKHSMDLSYALETVLNKQFPSSPQDVLGELQFAFVCFLLGNVYEAFEHWKRLLNLLCRSEAAMMKHHTLYINLISILYHQLGEIPADFFVDIVSQDNFLTSTLQVFFSSACSIAVDATLRKKAEKFQAHLTKKFRWDFAAEPEDCAPVVVELPEGIEMGXXMAENKGGGEAESGGGGSGSAPVTAGAAGPAAQEAEPPLTAVLVEEEEEEGGRAGAEGGAAGPDDGGVAAASSGSAQAASSPAASVGTGVAGGAVSTPAPAPASAPAPGPSAGPPPGPPASLLDTCAVCQQSLQSRREAEPKLLPCLHSFCLRCLPEPERQLSVPIPGGSNGDIQQVGVIRCPVCRQECRQIDLVDNYFVKDTSEAPSSSDEKSEQVCTSCEDNASAVGFCVECGEWLCKTCIEAHQRVKFTKDHLIRKKEDVSESVGASGQRPVFCPVHKQEQLKLFCETCDRLTCRDCQLLEHKEHRYQFLEEAFQNQKGAIENLLAKLLEKKNYVHFAATQVQNRIKEVNETNKRVEQEIKVAIFTLINEINKKGKSLLQQLENVTKERQMKLLQQQNDITGLSRQVKHVMNFTNWAIASGSSTALLYSKRLITFQLRHILKARCDPVPAANGAIRFHCDPTFWAKNVVNLGNLVIESKPAPGYTPNVVVGQVPPGTNHISKTPGQINLAQLRLQHMQQQVYAQKHQQLQQMRMQQPPAPVPTTTTTTQQHPRQAAPQMLQQQPPRLISVQTMQRGNMNCGAFQAHQMRLAQNAARIPGIPRHSGPQYSMMQPHLQRQHSNPGHAGPFPVVSVHNTTINPTSPTTATMANANRGPTSPSVTAIELIPSVTNPENLPSLPDIPPIQLEDAGSSSLDNLLSRYISGSHLPPQPTSTMNPSPGPSALSPGSSGLSNSHTPVRPPSTSSTGSRGSCGSSGRTAEKTSLSFKSDQVKVKQEPGTEDEICSFSGGVKQEKTEDGRRSACMLSSPESSLTPPLSTNLHLESELDALASLENHVKIEPADMNESCKQSGLSSLVNGKSPIRSLMHRSARIGGDGNNKDDDPNEDWCAVCQNGGDLLCCEKCPKVFHLTCHVPTLLSFPSGDWICTFCRDIGKPEVEYDCDNLQHSKKGKTAQGLSPVDQRKCERLLLYLYCHELSIEFQEPVPASIPNYYKIIKKPMDLSTVKKKLQKKHSQHYQIPDDFVADVRLIFKNCERFNEMMKVVQVYADTQEINLKADSEVAQAGKAVALYFEDKLTEIYSDRTFAPLPEFEQEEDDGEVTEDSDEDFIQPRRKRLKSDERPVHIKX"
## [6,] "MAAVQMDPELAKRLFFEGATVVILNMPKGTEFGIDYNSWEVGPKFRGVKMIPPGIHFLHYSSVDKANPKEVGPRMGFFLSLHQRGLTVLRWSTLREEVDLSPAPESEVEAMRANLQELDQFLGPYPYATLKKWISLTNFISEATVEKLQPENRQICAFSDVLPVLSMKHTKDRVGQNLPRCGIECKSYQEGLARLPEMKPRAGTEIRFSELPTQMFPEGATPAEITKHSMDLSYALETVLNKQFPSSPQDVLGELQFAFVCFLLGNVYEAFEHWKRLLNLLCRSEAAMMKHHTLYINLISILYHQLGEIPADFFVDIVSQDNFLTSTLQVFFSSACSIAVDATLRKKAEKFQAHLTKKFRWDFAAEPEDCAPVVVELPEGIEMGXXMAENKGGGEAESGGGGSGSAPVTAGAAGPAAQEAEPPLTAVLVEEEEEEGGRAGAEGGAAGPDDGGVAAASSGSAQAASSPAASVGTGVAGGAVSTPAPAPASAPAPGPSAGPPPGPPASLLDTCAVCQQSLQSRREAEPKLLPCLHSFCLRCLPEPERQLSVPIPGGSNGDIQQVGVIRCPVCRQECRQIDLVDNYFVKDTSEAPSSSDEKSEQVCTSCEDNASAVGFCVECGEWLCKTCIEAHQRVKFTKDHLIRKKEDVSESVGASGQRPVFCPVHKQEQLKLFCETCDRLTCRDCQLLEHKEHRYQFLEEAFQNQKGAIENLLAKLLEKKNYVHFAATQVQNRIKEVNETNKRVEQEIKVAIFTLINEINKKGKSLLQQLENVTKERQMKLLQQQNDITGLSRQVKHVMNFTNWAIASGSSTALLYSKRLITFQLRHILKARCDPVPAANGAIRFHCDPTFWAKNVVNLGNLVIESKPAPGYTPNVVVGQVPPGTNHISKTPGQINLAQLRLQHMQQQVYAQKHQQLQQMRMQQPPAPVPTTTTTTQQHPRQAAPQMLQQQPPRLISVQTMQRGNMNCGAFQAHQMRLAQNAARIPGIPRHSGPQYSMMQPHLQRQHSNPGHAGPFPVVSVHNTTINPTSPTTATMANANRGPTSPSVTAIELIPSVTNPENLPSLPDIPPIQLEDAGSSSLDNLLSRYISGSHLPPQPTSTMNPSPGPSALSPGSSGLSNSHTPVRPPSTSSTGSRGSCGSSGRTAEKTSLSFKSDQVKVKQEPGTEDEICSFSGGVKQEKTEDGRRSACMLSSPESSLTPPLSTNLHLESELDALASLENHVKIEPADMNESCKQSGLSSLVNGKSPIRSLMHRSARIGGDGNNKDDDPNEDWCAVCQNGGDLLCCEKCPKVFHLTCHVPTLLSFPSGDWICTFCRDIGKPEVEYDCDNLQHSKKGKTAQGLSPVDQRKCERLLLYLYCHELSIEFQEPVPASIPNYYKIIKKPMDLSTVKKKLQKKHSQHYQIPDDFVADVRLIFKNCERFNEMMKVVQVYADTQEINLKADSEVAQAGKAVALYFEDKLTEIYSDRTFAPLPEFEQEEDDGEVTEDSDEDFIQPRRKRLKSDERPVHIKX"
##      Mutant_Peptide                       Total_RNA Tumor_RNA_Ratio Tumor_RNA
## [1,] "DYNSWEVGPKFREQLKLFCETCDRLTCRDCQLLE" "3.03498" "NA"            "NA"     
## [2,] "DYNSWEVGPKFREQLKLFCETCDRLTCRDCQLLE" "3.03498" "NA"            "NA"     
## [3,] "DYNSWEVGPKFREQLKLFCETCDRLTCRDCQLLE" "3.03498" "NA"            "NA"     
## [4,] "DYNSWEVGPKFREQLKLFCETCDRLTCRDCQLLE" "3.03498" "NA"            "NA"     
## [5,] "DYNSWEVGPKFREQLKLFCETCDRLTCRDCQLLE" "3.03498" "NA"            "NA"     
## [6,] "DYNSWEVGPKFREQLKLFCETCDRLTCRDCQLLE" "3.03498" "NA"            "NA"     
##      Tumor_RNA_based_on_DNA MutRatio MutRatio_Min MutRatio_Max
## [1,] "NaN"                  "NA"     "NA"         "NA"        
## [2,] "NaN"                  "NA"     "NA"         "NA"        
## [3,] "NaN"                  "NA"     "NA"         "NA"        
## [4,] "NaN"                  "NA"     "NA"         "NA"        
## [5,] "NaN"                  "NA"     "NA"         "NA"        
## [6,] "NaN"                  "NA"     "NA"         "NA"
```

```
print(Export_Summary_IndelSV(Result_HLA1_SV, Mut_IC50_th = 500))
```

```
##              Num_All_Alteration        Num_Evaluated_Alteration 
##                              16                              16 
## Num_Alteration_Generating_NeoAg                 Num_All_Peptide 
##                              10                            3183 
##           Num_Evaluated_Peptide    Num_Peptide_Generating_NeoAg 
##                            3183                             119
```

```
print(Export_Summary_IndelSV_perFragments(Result_HLA1_SV, Mut_IC50_th = 500))
```

```
##                                                     DYNSWEVGPKFREQLKLFCETCDRLTCRDCQLLE-0_AAR2
## Num_Peptide_Per_Pep                                                                    57.000
## Num_Cond_Peptide_Per_Pep                                                               57.000
## Num_Rest_Peptide_Per_Pep                                                                0.000
## Num_Rest_Peptide_Per_Pep / Num_Cond_Peptide_Per_Pep                                     0.000
## -logP                                                                                   0.782
##                                                     GDPHELRNIFLQLSAVQEAQLKRLEVTRPRVLGSREQGQVPRMARQPPPPWVHAAFLLCLLSLGG-1_SLC25A12
## Num_Peptide_Per_Pep                                                                                                      237.000
## Num_Cond_Peptide_Per_Pep                                                                                                 237.000
## Num_Rest_Peptide_Per_Pep                                                                                                   6.000
## Num_Rest_Peptide_Per_Pep / Num_Cond_Peptide_Per_Pep                                                                        0.025
## -logP                                                                                                                      1.496
##                                                     TMAEKRQLFIEMLYRLX-2_DTNB
## Num_Peptide_Per_Pep                                                   24.000
## Num_Cond_Peptide_Per_Pep                                              24.000
## Num_Rest_Peptide_Per_Pep                                               6.000
## Num_Rest_Peptide_Per_Pep / Num_Cond_Peptide_Per_Pep                    0.250
## -logP                                                                  0.391
##                                                     TDWLSDCCFHPSERCRCTLYGHTDSVNSIEFFP-3_SPAG16
## Num_Peptide_Per_Pep                                                                     51.000
## Num_Cond_Peptide_Per_Pep                                                                51.000
## Num_Rest_Peptide_Per_Pep                                                                 0.000
## Num_Rest_Peptide_Per_Pep / Num_Cond_Peptide_Per_Pep                                      0.000
## -logP                                                                                    0.759
##                                                     ACLPGEEGTAERSCKRNLSTTTGSSERX-4_TMCC1
## Num_Peptide_Per_Pep                                                               90.000
## Num_Cond_Peptide_Per_Pep                                                          90.000
## Num_Rest_Peptide_Per_Pep                                                           0.000
## Num_Rest_Peptide_Per_Pep / Num_Cond_Peptide_Per_Pep                                0.000
## -logP                                                                              0.644
##                                                     MKSCKRNLSTTTGSSERX-5_TMCC1
## Num_Peptide_Per_Pep                                                     45.000
## Num_Cond_Peptide_Per_Pep                                                45.000
## Num_Rest_Peptide_Per_Pep                                                 0.000
## Num_Rest_Peptide_Per_Pep / Num_Cond_Peptide_Per_Pep                      0.000
## -logP                                                                    0.414
##                                                     VQRFSLRRQLSKSCKRNLSTTTGSSERX-6_TMCC1
## Num_Peptide_Per_Pep                                                               90.000
## Num_Cond_Peptide_Per_Pep                                                          90.000
## Num_Rest_Peptide_Per_Pep                                                           0.000
## Num_Rest_Peptide_Per_Pep / Num_Cond_Peptide_Per_Pep                                0.000
## -logP                                                                              0.644
##                                                     QRKGSYCWDFEESCKRNLSTTTGSSERX-7_TMCC1
## Num_Peptide_Per_Pep                                                               90.000
## Num_Cond_Peptide_Per_Pep                                                          90.000
## Num_Rest_Peptide_Per_Pep                                                           0.000
## Num_Rest_Peptide_Per_Pep / Num_Cond_Peptide_Per_Pep                                0.000
## -logP                                                                              0.644
##                                                     IRPKPFVIPGRSRTLECKQQKIILLRREMKETX-8_PDLIM3
## Num_Peptide_Per_Pep                                                                    120.000
## Num_Cond_Peptide_Per_Pep                                                               120.000
## Num_Rest_Peptide_Per_Pep                                                                 1.000
## Num_Rest_Peptide_Per_Pep / Num_Cond_Peptide_Per_Pep                                      0.008
## -logP                                                                                    0.759
##                                                     GQLRGLIPSSPQKTLECKQQKIILLRREMKETX-9_PDLIM3
## Num_Peptide_Per_Pep                                                                    120.000
## Num_Cond_Peptide_Per_Pep                                                               120.000
## Num_Rest_Peptide_Per_Pep                                                                 1.000
## Num_Rest_Peptide_Per_Pep / Num_Cond_Peptide_Per_Pep                                      0.008
## -logP                                                                                    0.759
##                                                     GDQYIVDMANTKRTLLWTRPMTVILGKPFCVMPKQQKTAHIGFLQHIPRLSPKPCLPKLNLMMRKQRMSQNGKNVKFEESHLRAVCMSGRGMGQVWVFFLCSX-10_WDR70
## Num_Peptide_Per_Pep                                                                                                                                          540.000
## Num_Cond_Peptide_Per_Pep                                                                                                                                     540.000
## Num_Rest_Peptide_Per_Pep                                                                                                                                      39.000
## Num_Rest_Peptide_Per_Pep / Num_Cond_Peptide_Per_Pep                                                                                                            0.072
## -logP                                                                                                                                                          2.370
##                                                     LCPASRALEEKKGHELLGGPSLGAPRPGSQERTGENTAACQDHRFWAGQTAGCGRERIPCRRRQSAYQVDGIGINFTQNLYPPEX-11_EGFR
## Num_Peptide_Per_Pep                                                                                                                       432.000
## Num_Cond_Peptide_Per_Pep                                                                                                                  432.000
## Num_Rest_Peptide_Per_Pep                                                                                                                    3.000
## Num_Rest_Peptide_Per_Pep / Num_Cond_Peptide_Per_Pep                                                                                         0.007
## -logP                                                                                                                                       1.956
##                                                     LTPMPEGLSQQQDMLKNTSKGHPDRLPLQMALT-12_ARHGEF10
## Num_Peptide_Per_Pep                                                                        51.000
## Num_Cond_Peptide_Per_Pep                                                                   51.000
## Num_Rest_Peptide_Per_Pep                                                                    0.000
## Num_Rest_Peptide_Per_Pep / Num_Cond_Peptide_Per_Pep                                         0.000
## -logP                                                                                       0.759
##                                                     PGVQGKKVRKPSPTTSASSSSWTSAGTSWPAGRRTGTATSGTATTTSVWPGCGTRMWSTQWSSVPRSRSCCSRPATTPPSKPGAPHAPCASSRHLAHGLAPSSPGLPARGAEVCWVHWSHRDPLRTSPGSVAFSRAGEVEMLIAVTPX-13_NOTCH1
## Num_Peptide_Per_Pep                                                                                                                                                                                        810.000
## Num_Cond_Peptide_Per_Pep                                                                                                                                                                                   810.000
## Num_Rest_Peptide_Per_Pep                                                                                                                                                                                    41.000
## Num_Rest_Peptide_Per_Pep / Num_Cond_Peptide_Per_Pep                                                                                                                                                          0.051
## -logP                                                                                                                                                                                                        3.405
##                                                     SRGRGNNNRKGREVTPL-14_UBAP2
## Num_Peptide_Per_Pep                                                     30.000
## Num_Cond_Peptide_Per_Pep                                                30.000
## Num_Rest_Peptide_Per_Pep                                                 0.000
## Num_Rest_Peptide_Per_Pep / Num_Cond_Peptide_Per_Pep                      0.000
## -logP                                                                    0.391
##                                                     IPSSQLAAKLLHMLTMRMLSKSATGRWX-15_RASGRP2
## Num_Peptide_Per_Pep                                                                  90.000
## Num_Cond_Peptide_Per_Pep                                                             90.000
## Num_Rest_Peptide_Per_Pep                                                             11.000
## Num_Rest_Peptide_Per_Pep / Num_Cond_Peptide_Per_Pep                                   0.122
## -logP                                                                                 0.644
##                                                     RKSPPEKKLRRYPPGQGATIGVDFMIKTVEINGE-16_GAB2
## Num_Peptide_Per_Pep                                                                     57.000
## Num_Cond_Peptide_Per_Pep                                                                57.000
## Num_Rest_Peptide_Per_Pep                                                                 0.000
## Num_Rest_Peptide_Per_Pep / Num_Cond_Peptide_Per_Pep                                      0.000
## -logP                                                                                    0.782
##                                                     MSVIVLPQPDEVLNLVQSYVTLRVPLYVSYVFH-17_MTRF1
## Num_Peptide_Per_Pep                                                                     51.000
## Num_Cond_Peptide_Per_Pep                                                                51.000
## Num_Rest_Peptide_Per_Pep                                                                 2.000
## Num_Rest_Peptide_Per_Pep / Num_Cond_Peptide_Per_Pep                                      0.039
## -logP                                                                                    0.759
##                                                     LWKCLRKPVLNDRNLQLHTDKGSFLKEKNKKLKKK-18_EVI2B
## Num_Peptide_Per_Pep                                                                       63.000
## Num_Cond_Peptide_Per_Pep                                                                  63.000
## Num_Rest_Peptide_Per_Pep                                                                   1.000
## Num_Rest_Peptide_Per_Pep / Num_Cond_Peptide_Per_Pep                                        0.016
## -logP                                                                                      0.805
##                                                     RTFLQSTDPGDTGGVFEGSTLCPSDPAFSRTDDPX-19_IKZF3
## Num_Peptide_Per_Pep                                                                      132.000
## Num_Cond_Peptide_Per_Pep                                                                 132.000
## Num_Rest_Peptide_Per_Pep                                                                   2.000
## Num_Rest_Peptide_Per_Pep / Num_Cond_Peptide_Per_Pep                                        0.015
## -logP                                                                                      0.805
##                                                     NVLMVHKRSHTGGVFEGSTLCPSDPAFSRTDDPX-20_IKZF3
## Num_Peptide_Per_Pep                                                                     126.000
## Num_Cond_Peptide_Per_Pep                                                                126.000
## Num_Rest_Peptide_Per_Pep                                                                  4.000
## Num_Rest_Peptide_Per_Pep / Num_Cond_Peptide_Per_Pep                                       0.032
## -logP                                                                                     0.782
##                                                     DALTGHLRTHSGGVFEGSTLCPSDPAFSRTDDPX-21_IKZF3
## Num_Peptide_Per_Pep                                                                     126.000
## Num_Cond_Peptide_Per_Pep                                                                126.000
## Num_Rest_Peptide_Per_Pep                                                                  4.000
## Num_Rest_Peptide_Per_Pep / Num_Cond_Peptide_Per_Pep                                       0.032
## -logP                                                                                     0.782
##                                                     FNVLMVHKRSHTGGVFEGSTLCPSDPAFSRTDDPX-22_IKZF3
## Num_Peptide_Per_Pep                                                                      132.000
## Num_Cond_Peptide_Per_Pep                                                                 132.000
## Num_Rest_Peptide_Per_Pep                                                                   4.000
## Num_Rest_Peptide_Per_Pep / Num_Cond_Peptide_Per_Pep                                        0.030
## -logP                                                                                      0.805
##                                                     EGPANEDEDIGGGVFEGSTLCPSDPAFSRTDDPX-23_IKZF3
## Num_Peptide_Per_Pep                                                                     126.000
## Num_Cond_Peptide_Per_Pep                                                                126.000
## Num_Rest_Peptide_Per_Pep                                                                  1.000
## Num_Rest_Peptide_Per_Pep / Num_Cond_Peptide_Per_Pep                                       0.008
## -logP                                                                                     0.782
##                                                     GEGPANEDEDIGGGVFEGSTLCPSDPAFSRTDDPX-24_IKZF3
## Num_Peptide_Per_Pep                                                                      132.000
## Num_Cond_Peptide_Per_Pep                                                                 132.000
## Num_Rest_Peptide_Per_Pep                                                                   1.000
## Num_Rest_Peptide_Per_Pep / Num_Cond_Peptide_Per_Pep                                        0.008
## -logP                                                                                      0.805
##                                                     RDALTGHLRTHSGGVFEGSTLCPSDPAFSRTDDPX-25_IKZF3
## Num_Peptide_Per_Pep                                                                      132.000
## Num_Cond_Peptide_Per_Pep                                                                 132.000
## Num_Rest_Peptide_Per_Pep                                                                   4.000
## Num_Rest_Peptide_Per_Pep / Num_Cond_Peptide_Per_Pep                                        0.030
## -logP                                                                                      0.805
```

```
print(head(Result_HLA2_SV))
```

```
#      HLA                     Pos Gene     Evaluated_Mutant_Peptide_Core
## [1,] "HLA-DPA10103-DPB10201" "6" "0_AAR2" "WEVGPKFRE"                  
## [2,] "HLA-DPA10103-DPB10201" "5" "0_AAR2" "WEVGPKFRE"                  
## [3,] "HLA-DPA10103-DPB10201" "4" "0_AAR2" "WEVGPKFRE"                  
## [4,] "HLA-DPA10103-DPB10201" "3" "0_AAR2" "WEVGPKFRE"                  
## [5,] "HLA-DPA10103-DPB10201" "2" "0_AAR2" "WEVGPKFRE"                  
## [6,] "HLA-DPA10103-DPB10201" "6" "0_AAR2" "KFREQLKLF"
##      Evaluated_Mutant_Peptide Mut_IC50  Mut_Rank Chr  NM_ID                
## [1,] "GIDYNSWEVGPKFRE"        "3510.19" "65.00"  "20" "NM_015511_NM_015906"
## [2,] "IDYNSWEVGPKFREQ"        "3548.49" "65.00"  "20" "NM_015511_NM_015906"
## [3,] "DYNSWEVGPKFREQL"        "3382.78" "65.00"  "20" "NM_015511_NM_015906"
## [4,] "YNSWEVGPKFREQLK"        "3259.68" "65.00"  "20" "NM_015511_NM_015906"
## [5,] "NSWEVGPKFREQLKL"        "3578.45" "65.00"  "20" "NM_015511_NM_015906"
## [6,] "SWEVGPKFREQLKLF"        "2000.52" "50.00"  "20" "NM_015511_NM_015906"
##      Change                     Ref Alt              Prob Mutation_Prob.
## [1,] "In_AAR2_exon_TRIM33_exon" "G" "G]1:115005805]" "0"  "0"           
## [2,] "In_AAR2_exon_TRIM33_exon" "G" "G]1:115005805]" "0"  "0"           
## [3,] "In_AAR2_exon_TRIM33_exon" "G" "G]1:115005805]" "0"  "0"           
## [4,] "In_AAR2_exon_TRIM33_exon" "G" "G]1:115005805]" "0"  "0"           
## [5,] "In_AAR2_exon_TRIM33_exon" "G" "G]1:115005805]" "0"  "0"           
## [6,] "In_AAR2_exon_TRIM33_exon" "G" "G]1:115005805]" "0"  "0"         
##      Exon_Start Exon_End   Mutation_Position Total_Depth Tumor_Depth
## [1,] "34824338" "34844863" "20_34827929"     "0"         "0"        
## [2,] "34824338" "34844863" "20_34827929"     "0"         "0"        
## [3,] "34824338" "34844863" "20_34827929"     "0"         "0"        
## [4,] "34824338" "34844863" "20_34827929"     "0"         "0"        
## [5,] "34824338" "34844863" "20_34827929"     "0"         "0"        
## [6,] "34824338" "34844863" "20_34827929"     "0"         "0"        
##      Wt_Peptide                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  
## [1,] "MAAVQMDPELAKRLFFEGATVVILNMPKGTEFGIDYNSWEVGPKFRGVKMIPPGIHFLHYSSVDKANPKEVGPRMGFFLSLHQRGLTVLRWSTLREEVDLSPAPESEVEAMRANLQELDQFLGPYPYATLKKWISLTNFISEATVEKLQPENRQICAFSDVLPVLSMKHTKDRVGQNLPRCGIECKSYQEGLARLPEMKPRAGTEIRFSELPTQMFPEGATPAEITKHSMDLSYALETVLNKQFPSSPQDVLGELQFAFVCFLLGNVYEAFEHWKRLLNLLCRSEAAMMKHHTLYINLISILYHQLGEIPADFFVDIVSQDNFLTSTLQVFFSSACSIAVDATLRKKAEKFQAHLTKKFRWDFAAEPEDCAPVVVELPEGIEMGXXMAENKGGGEAESGGGGSGSAPVTAGAAGPAAQEAEPPLTAVLVEEEEEEGGRAGAEGGAAGPDDGGVAAASSGSAQAASSPAASVGTGVAGGAVSTPAPAPASAPAPGPSAGPPPGPPASLLDTCAVCQQSLQSRREAEPKLLPCLHSFCLRCLPEPERQLSVPIPGGSNGDIQQVGVIRCPVCRQECRQIDLVDNYFVKDTSEAPSSSDEKSEQVCTSCEDNASAVGFCVECGEWLCKTCIEAHQRVKFTKDHLIRKKEDVSESVGASGQRPVFCPVHKQEQLKLFCETCDRLTCRDCQLLEHKEHRYQFLEEAFQNQKGAIENLLAKLLEKKNYVHFAATQVQNRIKEVNETNKRVEQEIKVAIFTLINEINKKGKSLLQQLENVTKERQMKLLQQQNDITGLSRQVKHVMNFTNWAIASGSSTALLYSKRLITFQLRHILKARCDPVPAANGAIRFHCDPTFWAKNVVNLGNLVIESKPAPGYTPNVVVGQVPPGTNHISKTPGQINLAQLRLQHMQQQVYAQKHQQLQQMRMQQPPAPVPTTTTTTQQHPRQAAPQMLQQQPPRLISVQTMQRGNMNCGAFQAHQMRLAQNAARIPGIPRHSGPQYSMMQPHLQRQHSNPGHAGPFPVVSVHNTTINPTSPTTATMANANRGPTSPSVTAIELIPSVTNPENLPSLPDIPPIQLEDAGSSSLDNLLSRYISGSHLPPQPTSTMNPSPGPSALSPGSSGLSNSHTPVRPPSTSSTGSRGSCGSSGRTAEKTSLSFKSDQVKVKQEPGTEDEICSFSGGVKQEKTEDGRRSACMLSSPESSLTPPLSTNLHLESELDALASLENHVKIEPADMNESCKQSGLSSLVNGKSPIRSLMHRSARIGGDGNNKDDDPNEDWCAVCQNGGDLLCCEKCPKVFHLTCHVPTLLSFPSGDWICTFCRDIGKPEVEYDCDNLQHSKKGKTAQGLSPVDQRKCERLLLYLYCHELSIEFQEPVPASIPNYYKIIKKPMDLSTVKKKLQKKHSQHYQIPDDFVADVRLIFKNCERFNEMMKVVQVYADTQEINLKADSEVAQAGKAVALYFEDKLTEIYSDRTFAPLPEFEQEEDDGEVTEDSDEDFIQPRRKRLKSDERPVHIKX"
## [2,] "MAAVQMDPELAKRLFFEGATVVILNMPKGTEFGIDYNSWEVGPKFRGVKMIPPGIHFLHYSSVDKANPKEVGPRMGFFLSLHQRGLTVLRWSTLREEVDLSPAPESEVEAMRANLQELDQFLGPYPYATLKKWISLTNFISEATVEKLQPENRQICAFSDVLPVLSMKHTKDRVGQNLPRCGIECKSYQEGLARLPEMKPRAGTEIRFSELPTQMFPEGATPAEITKHSMDLSYALETVLNKQFPSSPQDVLGELQFAFVCFLLGNVYEAFEHWKRLLNLLCRSEAAMMKHHTLYINLISILYHQLGEIPADFFVDIVSQDNFLTSTLQVFFSSACSIAVDATLRKKAEKFQAHLTKKFRWDFAAEPEDCAPVVVELPEGIEMGXXMAENKGGGEAESGGGGSGSAPVTAGAAGPAAQEAEPPLTAVLVEEEEEEGGRAGAEGGAAGPDDGGVAAASSGSAQAASSPAASVGTGVAGGAVSTPAPAPASAPAPGPSAGPPPGPPASLLDTCAVCQQSLQSRREAEPKLLPCLHSFCLRCLPEPERQLSVPIPGGSNGDIQQVGVIRCPVCRQECRQIDLVDNYFVKDTSEAPSSSDEKSEQVCTSCEDNASAVGFCVECGEWLCKTCIEAHQRVKFTKDHLIRKKEDVSESVGASGQRPVFCPVHKQEQLKLFCETCDRLTCRDCQLLEHKEHRYQFLEEAFQNQKGAIENLLAKLLEKKNYVHFAATQVQNRIKEVNETNKRVEQEIKVAIFTLINEINKKGKSLLQQLENVTKERQMKLLQQQNDITGLSRQVKHVMNFTNWAIASGSSTALLYSKRLITFQLRHILKARCDPVPAANGAIRFHCDPTFWAKNVVNLGNLVIESKPAPGYTPNVVVGQVPPGTNHISKTPGQINLAQLRLQHMQQQVYAQKHQQLQQMRMQQPPAPVPTTTTTTQQHPRQAAPQMLQQQPPRLISVQTMQRGNMNCGAFQAHQMRLAQNAARIPGIPRHSGPQYSMMQPHLQRQHSNPGHAGPFPVVSVHNTTINPTSPTTATMANANRGPTSPSVTAIELIPSVTNPENLPSLPDIPPIQLEDAGSSSLDNLLSRYISGSHLPPQPTSTMNPSPGPSALSPGSSGLSNSHTPVRPPSTSSTGSRGSCGSSGRTAEKTSLSFKSDQVKVKQEPGTEDEICSFSGGVKQEKTEDGRRSACMLSSPESSLTPPLSTNLHLESELDALASLENHVKIEPADMNESCKQSGLSSLVNGKSPIRSLMHRSARIGGDGNNKDDDPNEDWCAVCQNGGDLLCCEKCPKVFHLTCHVPTLLSFPSGDWICTFCRDIGKPEVEYDCDNLQHSKKGKTAQGLSPVDQRKCERLLLYLYCHELSIEFQEPVPASIPNYYKIIKKPMDLSTVKKKLQKKHSQHYQIPDDFVADVRLIFKNCERFNEMMKVVQVYADTQEINLKADSEVAQAGKAVALYFEDKLTEIYSDRTFAPLPEFEQEEDDGEVTEDSDEDFIQPRRKRLKSDERPVHIKX"
## [3,] "MAAVQMDPELAKRLFFEGATVVILNMPKGTEFGIDYNSWEVGPKFRGVKMIPPGIHFLHYSSVDKANPKEVGPRMGFFLSLHQRGLTVLRWSTLREEVDLSPAPESEVEAMRANLQELDQFLGPYPYATLKKWISLTNFISEATVEKLQPENRQICAFSDVLPVLSMKHTKDRVGQNLPRCGIECKSYQEGLARLPEMKPRAGTEIRFSELPTQMFPEGATPAEITKHSMDLSYALETVLNKQFPSSPQDVLGELQFAFVCFLLGNVYEAFEHWKRLLNLLCRSEAAMMKHHTLYINLISILYHQLGEIPADFFVDIVSQDNFLTSTLQVFFSSACSIAVDATLRKKAEKFQAHLTKKFRWDFAAEPEDCAPVVVELPEGIEMGXXMAENKGGGEAESGGGGSGSAPVTAGAAGPAAQEAEPPLTAVLVEEEEEEGGRAGAEGGAAGPDDGGVAAASSGSAQAASSPAASVGTGVAGGAVSTPAPAPASAPAPGPSAGPPPGPPASLLDTCAVCQQSLQSRREAEPKLLPCLHSFCLRCLPEPERQLSVPIPGGSNGDIQQVGVIRCPVCRQECRQIDLVDNYFVKDTSEAPSSSDEKSEQVCTSCEDNASAVGFCVECGEWLCKTCIEAHQRVKFTKDHLIRKKEDVSESVGASGQRPVFCPVHKQEQLKLFCETCDRLTCRDCQLLEHKEHRYQFLEEAFQNQKGAIENLLAKLLEKKNYVHFAATQVQNRIKEVNETNKRVEQEIKVAIFTLINEINKKGKSLLQQLENVTKERQMKLLQQQNDITGLSRQVKHVMNFTNWAIASGSSTALLYSKRLITFQLRHILKARCDPVPAANGAIRFHCDPTFWAKNVVNLGNLVIESKPAPGYTPNVVVGQVPPGTNHISKTPGQINLAQLRLQHMQQQVYAQKHQQLQQMRMQQPPAPVPTTTTTTQQHPRQAAPQMLQQQPPRLISVQTMQRGNMNCGAFQAHQMRLAQNAARIPGIPRHSGPQYSMMQPHLQRQHSNPGHAGPFPVVSVHNTTINPTSPTTATMANANRGPTSPSVTAIELIPSVTNPENLPSLPDIPPIQLEDAGSSSLDNLLSRYISGSHLPPQPTSTMNPSPGPSALSPGSSGLSNSHTPVRPPSTSSTGSRGSCGSSGRTAEKTSLSFKSDQVKVKQEPGTEDEICSFSGGVKQEKTEDGRRSACMLSSPESSLTPPLSTNLHLESELDALASLENHVKIEPADMNESCKQSGLSSLVNGKSPIRSLMHRSARIGGDGNNKDDDPNEDWCAVCQNGGDLLCCEKCPKVFHLTCHVPTLLSFPSGDWICTFCRDIGKPEVEYDCDNLQHSKKGKTAQGLSPVDQRKCERLLLYLYCHELSIEFQEPVPASIPNYYKIIKKPMDLSTVKKKLQKKHSQHYQIPDDFVADVRLIFKNCERFNEMMKVVQVYADTQEINLKADSEVAQAGKAVALYFEDKLTEIYSDRTFAPLPEFEQEEDDGEVTEDSDEDFIQPRRKRLKSDERPVHIKX"
## [4,] "MAAVQMDPELAKRLFFEGATVVILNMPKGTEFGIDYNSWEVGPKFRGVKMIPPGIHFLHYSSVDKANPKEVGPRMGFFLSLHQRGLTVLRWSTLREEVDLSPAPESEVEAMRANLQELDQFLGPYPYATLKKWISLTNFISEATVEKLQPENRQICAFSDVLPVLSMKHTKDRVGQNLPRCGIECKSYQEGLARLPEMKPRAGTEIRFSELPTQMFPEGATPAEITKHSMDLSYALETVLNKQFPSSPQDVLGELQFAFVCFLLGNVYEAFEHWKRLLNLLCRSEAAMMKHHTLYINLISILYHQLGEIPADFFVDIVSQDNFLTSTLQVFFSSACSIAVDATLRKKAEKFQAHLTKKFRWDFAAEPEDCAPVVVELPEGIEMGXXMAENKGGGEAESGGGGSGSAPVTAGAAGPAAQEAEPPLTAVLVEEEEEEGGRAGAEGGAAGPDDGGVAAASSGSAQAASSPAASVGTGVAGGAVSTPAPAPASAPAPGPSAGPPPGPPASLLDTCAVCQQSLQSRREAEPKLLPCLHSFCLRCLPEPERQLSVPIPGGSNGDIQQVGVIRCPVCRQECRQIDLVDNYFVKDTSEAPSSSDEKSEQVCTSCEDNASAVGFCVECGEWLCKTCIEAHQRVKFTKDHLIRKKEDVSESVGASGQRPVFCPVHKQEQLKLFCETCDRLTCRDCQLLEHKEHRYQFLEEAFQNQKGAIENLLAKLLEKKNYVHFAATQVQNRIKEVNETNKRVEQEIKVAIFTLINEINKKGKSLLQQLENVTKERQMKLLQQQNDITGLSRQVKHVMNFTNWAIASGSSTALLYSKRLITFQLRHILKARCDPVPAANGAIRFHCDPTFWAKNVVNLGNLVIESKPAPGYTPNVVVGQVPPGTNHISKTPGQINLAQLRLQHMQQQVYAQKHQQLQQMRMQQPPAPVPTTTTTTQQHPRQAAPQMLQQQPPRLISVQTMQRGNMNCGAFQAHQMRLAQNAARIPGIPRHSGPQYSMMQPHLQRQHSNPGHAGPFPVVSVHNTTINPTSPTTATMANANRGPTSPSVTAIELIPSVTNPENLPSLPDIPPIQLEDAGSSSLDNLLSRYISGSHLPPQPTSTMNPSPGPSALSPGSSGLSNSHTPVRPPSTSSTGSRGSCGSSGRTAEKTSLSFKSDQVKVKQEPGTEDEICSFSGGVKQEKTEDGRRSACMLSSPESSLTPPLSTNLHLESELDALASLENHVKIEPADMNESCKQSGLSSLVNGKSPIRSLMHRSARIGGDGNNKDDDPNEDWCAVCQNGGDLLCCEKCPKVFHLTCHVPTLLSFPSGDWICTFCRDIGKPEVEYDCDNLQHSKKGKTAQGLSPVDQRKCERLLLYLYCHELSIEFQEPVPASIPNYYKIIKKPMDLSTVKKKLQKKHSQHYQIPDDFVADVRLIFKNCERFNEMMKVVQVYADTQEINLKADSEVAQAGKAVALYFEDKLTEIYSDRTFAPLPEFEQEEDDGEVTEDSDEDFIQPRRKRLKSDERPVHIKX"
## [5,] "MAAVQMDPELAKRLFFEGATVVILNMPKGTEFGIDYNSWEVGPKFRGVKMIPPGIHFLHYSSVDKANPKEVGPRMGFFLSLHQRGLTVLRWSTLREEVDLSPAPESEVEAMRANLQELDQFLGPYPYATLKKWISLTNFISEATVEKLQPENRQICAFSDVLPVLSMKHTKDRVGQNLPRCGIECKSYQEGLARLPEMKPRAGTEIRFSELPTQMFPEGATPAEITKHSMDLSYALETVLNKQFPSSPQDVLGELQFAFVCFLLGNVYEAFEHWKRLLNLLCRSEAAMMKHHTLYINLISILYHQLGEIPADFFVDIVSQDNFLTSTLQVFFSSACSIAVDATLRKKAEKFQAHLTKKFRWDFAAEPEDCAPVVVELPEGIEMGXXMAENKGGGEAESGGGGSGSAPVTAGAAGPAAQEAEPPLTAVLVEEEEEEGGRAGAEGGAAGPDDGGVAAASSGSAQAASSPAASVGTGVAGGAVSTPAPAPASAPAPGPSAGPPPGPPASLLDTCAVCQQSLQSRREAEPKLLPCLHSFCLRCLPEPERQLSVPIPGGSNGDIQQVGVIRCPVCRQECRQIDLVDNYFVKDTSEAPSSSDEKSEQVCTSCEDNASAVGFCVECGEWLCKTCIEAHQRVKFTKDHLIRKKEDVSESVGASGQRPVFCPVHKQEQLKLFCETCDRLTCRDCQLLEHKEHRYQFLEEAFQNQKGAIENLLAKLLEKKNYVHFAATQVQNRIKEVNETNKRVEQEIKVAIFTLINEINKKGKSLLQQLENVTKERQMKLLQQQNDITGLSRQVKHVMNFTNWAIASGSSTALLYSKRLITFQLRHILKARCDPVPAANGAIRFHCDPTFWAKNVVNLGNLVIESKPAPGYTPNVVVGQVPPGTNHISKTPGQINLAQLRLQHMQQQVYAQKHQQLQQMRMQQPPAPVPTTTTTTQQHPRQAAPQMLQQQPPRLISVQTMQRGNMNCGAFQAHQMRLAQNAARIPGIPRHSGPQYSMMQPHLQRQHSNPGHAGPFPVVSVHNTTINPTSPTTATMANANRGPTSPSVTAIELIPSVTNPENLPSLPDIPPIQLEDAGSSSLDNLLSRYISGSHLPPQPTSTMNPSPGPSALSPGSSGLSNSHTPVRPPSTSSTGSRGSCGSSGRTAEKTSLSFKSDQVKVKQEPGTEDEICSFSGGVKQEKTEDGRRSACMLSSPESSLTPPLSTNLHLESELDALASLENHVKIEPADMNESCKQSGLSSLVNGKSPIRSLMHRSARIGGDGNNKDDDPNEDWCAVCQNGGDLLCCEKCPKVFHLTCHVPTLLSFPSGDWICTFCRDIGKPEVEYDCDNLQHSKKGKTAQGLSPVDQRKCERLLLYLYCHELSIEFQEPVPASIPNYYKIIKKPMDLSTVKKKLQKKHSQHYQIPDDFVADVRLIFKNCERFNEMMKVVQVYADTQEINLKADSEVAQAGKAVALYFEDKLTEIYSDRTFAPLPEFEQEEDDGEVTEDSDEDFIQPRRKRLKSDERPVHIKX"
## [6,] "MAAVQMDPELAKRLFFEGATVVILNMPKGTEFGIDYNSWEVGPKFRGVKMIPPGIHFLHYSSVDKANPKEVGPRMGFFLSLHQRGLTVLRWSTLREEVDLSPAPESEVEAMRANLQELDQFLGPYPYATLKKWISLTNFISEATVEKLQPENRQICAFSDVLPVLSMKHTKDRVGQNLPRCGIECKSYQEGLARLPEMKPRAGTEIRFSELPTQMFPEGATPAEITKHSMDLSYALETVLNKQFPSSPQDVLGELQFAFVCFLLGNVYEAFEHWKRLLNLLCRSEAAMMKHHTLYINLISILYHQLGEIPADFFVDIVSQDNFLTSTLQVFFSSACSIAVDATLRKKAEKFQAHLTKKFRWDFAAEPEDCAPVVVELPEGIEMGXXMAENKGGGEAESGGGGSGSAPVTAGAAGPAAQEAEPPLTAVLVEEEEEEGGRAGAEGGAAGPDDGGVAAASSGSAQAASSPAASVGTGVAGGAVSTPAPAPASAPAPGPSAGPPPGPPASLLDTCAVCQQSLQSRREAEPKLLPCLHSFCLRCLPEPERQLSVPIPGGSNGDIQQVGVIRCPVCRQECRQIDLVDNYFVKDTSEAPSSSDEKSEQVCTSCEDNASAVGFCVECGEWLCKTCIEAHQRVKFTKDHLIRKKEDVSESVGASGQRPVFCPVHKQEQLKLFCETCDRLTCRDCQLLEHKEHRYQFLEEAFQNQKGAIENLLAKLLEKKNYVHFAATQVQNRIKEVNETNKRVEQEIKVAIFTLINEINKKGKSLLQQLENVTKERQMKLLQQQNDITGLSRQVKHVMNFTNWAIASGSSTALLYSKRLITFQLRHILKARCDPVPAANGAIRFHCDPTFWAKNVVNLGNLVIESKPAPGYTPNVVVGQVPPGTNHISKTPGQINLAQLRLQHMQQQVYAQKHQQLQQMRMQQPPAPVPTTTTTTQQHPRQAAPQMLQQQPPRLISVQTMQRGNMNCGAFQAHQMRLAQNAARIPGIPRHSGPQYSMMQPHLQRQHSNPGHAGPFPVVSVHNTTINPTSPTTATMANANRGPTSPSVTAIELIPSVTNPENLPSLPDIPPIQLEDAGSSSLDNLLSRYISGSHLPPQPTSTMNPSPGPSALSPGSSGLSNSHTPVRPPSTSSTGSRGSCGSSGRTAEKTSLSFKSDQVKVKQEPGTEDEICSFSGGVKQEKTEDGRRSACMLSSPESSLTPPLSTNLHLESELDALASLENHVKIEPADMNESCKQSGLSSLVNGKSPIRSLMHRSARIGGDGNNKDDDPNEDWCAVCQNGGDLLCCEKCPKVFHLTCHVPTLLSFPSGDWICTFCRDIGKPEVEYDCDNLQHSKKGKTAQGLSPVDQRKCERLLLYLYCHELSIEFQEPVPASIPNYYKIIKKPMDLSTVKKKLQKKHSQHYQIPDDFVADVRLIFKNCERFNEMMKVVQVYADTQEINLKADSEVAQAGKAVALYFEDKLTEIYSDRTFAPLPEFEQEEDDGEVTEDSDEDFIQPRRKRLKSDERPVHIKX"
##      Mutant_Peptide                           Total_RNA Tumor_RNA_Ratio
## [1,] "GIDYNSWEVGPKFREQLKLFCETCDRLTCRDCQLLEHK" "3.03498" "NA"           
## [2,] "GIDYNSWEVGPKFREQLKLFCETCDRLTCRDCQLLEHK" "3.03498" "NA"           
## [3,] "GIDYNSWEVGPKFREQLKLFCETCDRLTCRDCQLLEHK" "3.03498" "NA"           
## [4,] "GIDYNSWEVGPKFREQLKLFCETCDRLTCRDCQLLEHK" "3.03498" "NA"           
## [5,] "GIDYNSWEVGPKFREQLKLFCETCDRLTCRDCQLLEHK" "3.03498" "NA"           
## [6,] "GIDYNSWEVGPKFREQLKLFCETCDRLTCRDCQLLEHK" "3.03498" "NA"           
##      Tumor_RNA Tumor_RNA_based_on_DNA MutRatio MutRatio_Min MutRatio_Max
## [1,] "NA"      "NaN"                  "NA"     "NA"         "NA"        
## [2,] "NA"      "NaN"                  "NA"     "NA"         "NA"        
## [3,] "NA"      "NaN"                  "NA"     "NA"         "NA"        
## [4,] "NA"      "NaN"                  "NA"     "NA"         "NA"        
## [5,] "NA"      "NaN"                  "NA"     "NA"         "NA"        
## [6,] "NA"      "NaN"                  "NA"     "NA"         "NA"
```

```
print(Export_Summary_IndelSV(Result_HLA2_SV, Mut_IC50_th = 500))
```

```
##              Num_All_Alteration        Num_Evaluated_Alteration 
##                              16                              16 
## Num_Alteration_Generating_NeoAg                 Num_All_Peptide 
##                              12                             588 
##           Num_Evaluated_Peptide    Num_Peptide_Generating_NeoAg 
##                             588                             298
```

  print(head(Result_HLA1_Seq))
```
##      HLA           Pos Gene           Evaluated_Mutant_Peptide_Core
## [1,] "HLA-A*02:01" "1" "0_atggcagaag" "MAEDDPYLGRPEK"              
## [2,] "HLA-A*02:01" "2" "0_atggcagaag" "AEDDPYLGRPEKM"              
## [3,] "HLA-A*02:01" "3" "0_atggcagaag" "YLGRPEKMF"                  
## [4,] "HLA-A*02:01" "4" "0_atggcagaag" "YLGRPEKMFH"                 
## [5,] "HLA-A*02:01" "5" "0_atggcagaag" "YLGRPEKMFHL"                
## [6,] "HLA-A*02:01" "6" "0_atggcagaag" "YLGRPEKMFHL"                
##      Evaluated_Mutant_Peptide Mut_IC50  Mut_Rank  Chr NM_ID ReadingFrame
## [1,] "MAEDDPYLGRPEK"          "34529.5" "52.6059" "0" ""    "1"         
## [2,] "AEDDPYLGRPEKM"          "37395.3" "61.4560" "0" ""    "1"         
## [3,] "EDDPYLGRPEKMF"          "41328.0" "76.2274" "0" ""    "1"         
## [4,] "DDPYLGRPEKMFH"          "42864.6" "82.5363" "0" ""    "1"         
## [5,] "DPYLGRPEKMFHL"          "3345.9"  "7.5921"  "0" ""    "1"         
## [6,] "PYLGRPEKMFHLD"          "19715.7" "24.8190" "0" ""    "1"
##      SequenceNumber Chrs        NM_IDs                   GeneIDs      
## [1,] "1"            "chr4;chr4" "NM_001165412;NM_003998" "NFKB1;NFKB1"
## [2,] "1"            "chr4;chr4" "NM_001165412;NM_003998" "NFKB1;NFKB1"
## [3,] "1"            "chr4;chr4" "NM_001165412;NM_003998" "NFKB1;NFKB1"
## [4,] "1"            "chr4;chr4" "NM_001165412;NM_003998" "NFKB1;NFKB1"
## [5,] "1"            "chr4;chr4" "NM_001165412;NM_003998" "NFKB1;NFKB1"
## [6,] "1"            "chr4;chr4" "NM_001165412;NM_003998" "NFKB1;NFKB1"
##      Exon_Starts           Exon_Ends             GroupID NumOfPeptides
## [1,] "103422515;103422515" "103538459;103538459" "0_1"   "27"         
## [2,] "103422515;103422515" "103538459;103538459" "0_1"   "27"         
## [3,] "103422515;103422515" "103538459;103538459" "0_1"   "27"         
## [4,] "103422515;103422515" "103538459;103538459" "0_1"   "27"         
## [5,] "103422515;103422515" "103538459;103538459" "0_1"   "27"         
## [6,] "103422515;103422515" "103538459;103538459" "0_1"   "27"         
##      NumOfStops
## [1,] "0"       
## [2,] "0"       
## [3,] "0"       
## [4,] "0"       
## [5,] "0"       
## [6,] "0"       
##      Wt_Peptide                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            
## [1,] "MAEDDPYLGRPEQMFHLDPSLTHTIFNPEVFQPQMALPTDGPYLQILEQPKQRGFRFRYVCEGPSHGGLPGASSEKNKKSYPQVKICNYVGPAKVIVQLVTNGKNIHLHAHSLVGKHCEDGICTVTAGPKDMVVGFANLGILHVTKKKVFETLEARMTEACIRGYNPGLLVHPDLAYLQAEGGGDRQLGDREKELIRQAALQQTKEMDLSVVRLMFTAFLPDSTGSFTRRLEPVVSDAIYDSKAPNASNLKIVRMDRTAGCVTGGEEIYLLCDKVQKDDIQIRFYEEEENGGVWEGFGDFSPTDVHRQFAIVFKTPKYKDINITKPASVFVQLRRKSDLETSEPKPFLYYPEIKDKEEVQRKRQKLMPNFSDSFGGGSGAGAGGGGMFGSGGGGGGTGSTGPGYSFPHYGFPTYGGITFHPGTTKSNAGMKHGTMDTESKKDPEGCDKSDDKNTVNLFGKVIETTEQDQEPSEATVGNGEVTLTYATGTKEESAGVQDNLFLEKAMQLAKRHANALFDYAVTGDVKMLLAVQRHLTAVQDENGDSVLHLAIIHLHSQLVRDLLEVTSGLISDDIINMRNDLYQTPLHLAVITKQEDVVEDLLRAGADLSLLDRLGNSVLHLAAKEGHDKVLSILLKHKKAALLLDHPNGDGLNAIHLAMMSNSLPCLLLLVAAGADVNAQEQKSGRTALHLAVEHDNISLAGCLLLEGDAHVDSTTYDGTTPLHIAAGRGSTRLAALLKAAGADPLVENFEPLYDLDDSWENAGEDEGVVPGTTPLDMATSWQVFDILNGKPYEPEFTSDDLLAQGDMKQLAEDVKLQLYKLLEIPDPDKNWATLAQKLGLGILNNAFRLSPAPSKTLMDNYEVSGGTVRELVEALRQMGYTEAIEVIQAASSPVKTTSQAHSLPLSPASTRQQIDELRDSDSVCDSGVETSFRKLSFTESLTSGASLLTLNKMPHDYGQEGPLEGKIXXMAEDDPYLGRPEQMFHLDPSLTHTIFNPEVFQPQMALPTADGPYLQILEQPKQRGFRFRYVCEGPSHGGLPGASSEKNKKSYPQVKICNYVGPAKVIVQLVTNGKNIHLHAHSLVGKHCEDGICTVTAGPKDMVVGFANLGILHVTKKKVFETLEARMTEACIRGYNPGLLVHPDLAYLQAEGGGDRQLGDREKELIRQAALQQTKEMDLSVVRLMFTAFLPDSTGSFTRRLEPVVSDAIYDSKAPNASNLKIVRMDRTAGCVTGGEEIYLLCDKVQKDDIQIRFYEEEENGGVWEGFGDFSPTDVHRQFAIVFKTPKYKDINITKPASVFVQLRRKSDLETSEPKPFLYYPEIKDKEEVQRKRQKLMPNFSDSFGGGSGAGAGGGGMFGSGGGGGGTGSTGPGYSFPHYGFPTYGGITFHPGTTKSNAGMKHGTMDTESKKDPEGCDKSDDKNTVNLFGKVIETTEQDQEPSEATVGNGEVTLTYATGTKEESAGVQDNLFLEKAMQLAKRHANALFDYAVTGDVKMLLAVQRHLTAVQDENGDSVLHLAIIHLHSQLVRDLLEVTSGLISDDIINMRNDLYQTPLHLAVITKQEDVVEDLLRAGADLSLLDRLGNSVLHLAAKEGHDKVLSILLKHKKAALLLDHPNGDGLNAIHLAMMSNSLPCLLLLVAAGADVNAQEQKSGRTALHLAVEHDNISLAGCLLLEGDAHVDSTTYDGTTPLHIAAGRGSTRLAALLKAAGADPLVENFEPLYDLDDSWENAGEDEGVVPGTTPLDMATSWQVFDILNGKPYEPEFTSDDLLAQGDMKQLAEDVKLQLYKLLEIPDPDKNWATLAQKLGLGILNNAFRLSPAPSKTLMDNYEVSGGTVRELVEALRQMGYTEAIEVIQAASSPVKTTSQAHSLPLSPASTRQQIDELRDSDSVCDSGVETSFRKLSFTESLTSGASLLTLNKMPHDYGQEGPLEGKIX"
## [2,] "MAEDDPYLGRPEQMFHLDPSLTHTIFNPEVFQPQMALPTDGPYLQILEQPKQRGFRFRYVCEGPSHGGLPGASSEKNKKSYPQVKICNYVGPAKVIVQLVTNGKNIHLHAHSLVGKHCEDGICTVTAGPKDMVVGFANLGILHVTKKKVFETLEARMTEACIRGYNPGLLVHPDLAYLQAEGGGDRQLGDREKELIRQAALQQTKEMDLSVVRLMFTAFLPDSTGSFTRRLEPVVSDAIYDSKAPNASNLKIVRMDRTAGCVTGGEEIYLLCDKVQKDDIQIRFYEEEENGGVWEGFGDFSPTDVHRQFAIVFKTPKYKDINITKPASVFVQLRRKSDLETSEPKPFLYYPEIKDKEEVQRKRQKLMPNFSDSFGGGSGAGAGGGGMFGSGGGGGGTGSTGPGYSFPHYGFPTYGGITFHPGTTKSNAGMKHGTMDTESKKDPEGCDKSDDKNTVNLFGKVIETTEQDQEPSEATVGNGEVTLTYATGTKEESAGVQDNLFLEKAMQLAKRHANALFDYAVTGDVKMLLAVQRHLTAVQDENGDSVLHLAIIHLHSQLVRDLLEVTSGLISDDIINMRNDLYQTPLHLAVITKQEDVVEDLLRAGADLSLLDRLGNSVLHLAAKEGHDKVLSILLKHKKAALLLDHPNGDGLNAIHLAMMSNSLPCLLLLVAAGADVNAQEQKSGRTALHLAVEHDNISLAGCLLLEGDAHVDSTTYDGTTPLHIAAGRGSTRLAALLKAAGADPLVENFEPLYDLDDSWENAGEDEGVVPGTTPLDMATSWQVFDILNGKPYEPEFTSDDLLAQGDMKQLAEDVKLQLYKLLEIPDPDKNWATLAQKLGLGILNNAFRLSPAPSKTLMDNYEVSGGTVRELVEALRQMGYTEAIEVIQAASSPVKTTSQAHSLPLSPASTRQQIDELRDSDSVCDSGVETSFRKLSFTESLTSGASLLTLNKMPHDYGQEGPLEGKIXXMAEDDPYLGRPEQMFHLDPSLTHTIFNPEVFQPQMALPTADGPYLQILEQPKQRGFRFRYVCEGPSHGGLPGASSEKNKKSYPQVKICNYVGPAKVIVQLVTNGKNIHLHAHSLVGKHCEDGICTVTAGPKDMVVGFANLGILHVTKKKVFETLEARMTEACIRGYNPGLLVHPDLAYLQAEGGGDRQLGDREKELIRQAALQQTKEMDLSVVRLMFTAFLPDSTGSFTRRLEPVVSDAIYDSKAPNASNLKIVRMDRTAGCVTGGEEIYLLCDKVQKDDIQIRFYEEEENGGVWEGFGDFSPTDVHRQFAIVFKTPKYKDINITKPASVFVQLRRKSDLETSEPKPFLYYPEIKDKEEVQRKRQKLMPNFSDSFGGGSGAGAGGGGMFGSGGGGGGTGSTGPGYSFPHYGFPTYGGITFHPGTTKSNAGMKHGTMDTESKKDPEGCDKSDDKNTVNLFGKVIETTEQDQEPSEATVGNGEVTLTYATGTKEESAGVQDNLFLEKAMQLAKRHANALFDYAVTGDVKMLLAVQRHLTAVQDENGDSVLHLAIIHLHSQLVRDLLEVTSGLISDDIINMRNDLYQTPLHLAVITKQEDVVEDLLRAGADLSLLDRLGNSVLHLAAKEGHDKVLSILLKHKKAALLLDHPNGDGLNAIHLAMMSNSLPCLLLLVAAGADVNAQEQKSGRTALHLAVEHDNISLAGCLLLEGDAHVDSTTYDGTTPLHIAAGRGSTRLAALLKAAGADPLVENFEPLYDLDDSWENAGEDEGVVPGTTPLDMATSWQVFDILNGKPYEPEFTSDDLLAQGDMKQLAEDVKLQLYKLLEIPDPDKNWATLAQKLGLGILNNAFRLSPAPSKTLMDNYEVSGGTVRELVEALRQMGYTEAIEVIQAASSPVKTTSQAHSLPLSPASTRQQIDELRDSDSVCDSGVETSFRKLSFTESLTSGASLLTLNKMPHDYGQEGPLEGKIX"
## [3,] "MAEDDPYLGRPEQMFHLDPSLTHTIFNPEVFQPQMALPTDGPYLQILEQPKQRGFRFRYVCEGPSHGGLPGASSEKNKKSYPQVKICNYVGPAKVIVQLVTNGKNIHLHAHSLVGKHCEDGICTVTAGPKDMVVGFANLGILHVTKKKVFETLEARMTEACIRGYNPGLLVHPDLAYLQAEGGGDRQLGDREKELIRQAALQQTKEMDLSVVRLMFTAFLPDSTGSFTRRLEPVVSDAIYDSKAPNASNLKIVRMDRTAGCVTGGEEIYLLCDKVQKDDIQIRFYEEEENGGVWEGFGDFSPTDVHRQFAIVFKTPKYKDINITKPASVFVQLRRKSDLETSEPKPFLYYPEIKDKEEVQRKRQKLMPNFSDSFGGGSGAGAGGGGMFGSGGGGGGTGSTGPGYSFPHYGFPTYGGITFHPGTTKSNAGMKHGTMDTESKKDPEGCDKSDDKNTVNLFGKVIETTEQDQEPSEATVGNGEVTLTYATGTKEESAGVQDNLFLEKAMQLAKRHANALFDYAVTGDVKMLLAVQRHLTAVQDENGDSVLHLAIIHLHSQLVRDLLEVTSGLISDDIINMRNDLYQTPLHLAVITKQEDVVEDLLRAGADLSLLDRLGNSVLHLAAKEGHDKVLSILLKHKKAALLLDHPNGDGLNAIHLAMMSNSLPCLLLLVAAGADVNAQEQKSGRTALHLAVEHDNISLAGCLLLEGDAHVDSTTYDGTTPLHIAAGRGSTRLAALLKAAGADPLVENFEPLYDLDDSWENAGEDEGVVPGTTPLDMATSWQVFDILNGKPYEPEFTSDDLLAQGDMKQLAEDVKLQLYKLLEIPDPDKNWATLAQKLGLGILNNAFRLSPAPSKTLMDNYEVSGGTVRELVEALRQMGYTEAIEVIQAASSPVKTTSQAHSLPLSPASTRQQIDELRDSDSVCDSGVETSFRKLSFTESLTSGASLLTLNKMPHDYGQEGPLEGKIXXMAEDDPYLGRPEQMFHLDPSLTHTIFNPEVFQPQMALPTADGPYLQILEQPKQRGFRFRYVCEGPSHGGLPGASSEKNKKSYPQVKICNYVGPAKVIVQLVTNGKNIHLHAHSLVGKHCEDGICTVTAGPKDMVVGFANLGILHVTKKKVFETLEARMTEACIRGYNPGLLVHPDLAYLQAEGGGDRQLGDREKELIRQAALQQTKEMDLSVVRLMFTAFLPDSTGSFTRRLEPVVSDAIYDSKAPNASNLKIVRMDRTAGCVTGGEEIYLLCDKVQKDDIQIRFYEEEENGGVWEGFGDFSPTDVHRQFAIVFKTPKYKDINITKPASVFVQLRRKSDLETSEPKPFLYYPEIKDKEEVQRKRQKLMPNFSDSFGGGSGAGAGGGGMFGSGGGGGGTGSTGPGYSFPHYGFPTYGGITFHPGTTKSNAGMKHGTMDTESKKDPEGCDKSDDKNTVNLFGKVIETTEQDQEPSEATVGNGEVTLTYATGTKEESAGVQDNLFLEKAMQLAKRHANALFDYAVTGDVKMLLAVQRHLTAVQDENGDSVLHLAIIHLHSQLVRDLLEVTSGLISDDIINMRNDLYQTPLHLAVITKQEDVVEDLLRAGADLSLLDRLGNSVLHLAAKEGHDKVLSILLKHKKAALLLDHPNGDGLNAIHLAMMSNSLPCLLLLVAAGADVNAQEQKSGRTALHLAVEHDNISLAGCLLLEGDAHVDSTTYDGTTPLHIAAGRGSTRLAALLKAAGADPLVENFEPLYDLDDSWENAGEDEGVVPGTTPLDMATSWQVFDILNGKPYEPEFTSDDLLAQGDMKQLAEDVKLQLYKLLEIPDPDKNWATLAQKLGLGILNNAFRLSPAPSKTLMDNYEVSGGTVRELVEALRQMGYTEAIEVIQAASSPVKTTSQAHSLPLSPASTRQQIDELRDSDSVCDSGVETSFRKLSFTESLTSGASLLTLNKMPHDYGQEGPLEGKIX"
## [4,] "MAEDDPYLGRPEQMFHLDPSLTHTIFNPEVFQPQMALPTDGPYLQILEQPKQRGFRFRYVCEGPSHGGLPGASSEKNKKSYPQVKICNYVGPAKVIVQLVTNGKNIHLHAHSLVGKHCEDGICTVTAGPKDMVVGFANLGILHVTKKKVFETLEARMTEACIRGYNPGLLVHPDLAYLQAEGGGDRQLGDREKELIRQAALQQTKEMDLSVVRLMFTAFLPDSTGSFTRRLEPVVSDAIYDSKAPNASNLKIVRMDRTAGCVTGGEEIYLLCDKVQKDDIQIRFYEEEENGGVWEGFGDFSPTDVHRQFAIVFKTPKYKDINITKPASVFVQLRRKSDLETSEPKPFLYYPEIKDKEEVQRKRQKLMPNFSDSFGGGSGAGAGGGGMFGSGGGGGGTGSTGPGYSFPHYGFPTYGGITFHPGTTKSNAGMKHGTMDTESKKDPEGCDKSDDKNTVNLFGKVIETTEQDQEPSEATVGNGEVTLTYATGTKEESAGVQDNLFLEKAMQLAKRHANALFDYAVTGDVKMLLAVQRHLTAVQDENGDSVLHLAIIHLHSQLVRDLLEVTSGLISDDIINMRNDLYQTPLHLAVITKQEDVVEDLLRAGADLSLLDRLGNSVLHLAAKEGHDKVLSILLKHKKAALLLDHPNGDGLNAIHLAMMSNSLPCLLLLVAAGADVNAQEQKSGRTALHLAVEHDNISLAGCLLLEGDAHVDSTTYDGTTPLHIAAGRGSTRLAALLKAAGADPLVENFEPLYDLDDSWENAGEDEGVVPGTTPLDMATSWQVFDILNGKPYEPEFTSDDLLAQGDMKQLAEDVKLQLYKLLEIPDPDKNWATLAQKLGLGILNNAFRLSPAPSKTLMDNYEVSGGTVRELVEALRQMGYTEAIEVIQAASSPVKTTSQAHSLPLSPASTRQQIDELRDSDSVCDSGVETSFRKLSFTESLTSGASLLTLNKMPHDYGQEGPLEGKIXXMAEDDPYLGRPEQMFHLDPSLTHTIFNPEVFQPQMALPTADGPYLQILEQPKQRGFRFRYVCEGPSHGGLPGASSEKNKKSYPQVKICNYVGPAKVIVQLVTNGKNIHLHAHSLVGKHCEDGICTVTAGPKDMVVGFANLGILHVTKKKVFETLEARMTEACIRGYNPGLLVHPDLAYLQAEGGGDRQLGDREKELIRQAALQQTKEMDLSVVRLMFTAFLPDSTGSFTRRLEPVVSDAIYDSKAPNASNLKIVRMDRTAGCVTGGEEIYLLCDKVQKDDIQIRFYEEEENGGVWEGFGDFSPTDVHRQFAIVFKTPKYKDINITKPASVFVQLRRKSDLETSEPKPFLYYPEIKDKEEVQRKRQKLMPNFSDSFGGGSGAGAGGGGMFGSGGGGGGTGSTGPGYSFPHYGFPTYGGITFHPGTTKSNAGMKHGTMDTESKKDPEGCDKSDDKNTVNLFGKVIETTEQDQEPSEATVGNGEVTLTYATGTKEESAGVQDNLFLEKAMQLAKRHANALFDYAVTGDVKMLLAVQRHLTAVQDENGDSVLHLAIIHLHSQLVRDLLEVTSGLISDDIINMRNDLYQTPLHLAVITKQEDVVEDLLRAGADLSLLDRLGNSVLHLAAKEGHDKVLSILLKHKKAALLLDHPNGDGLNAIHLAMMSNSLPCLLLLVAAGADVNAQEQKSGRTALHLAVEHDNISLAGCLLLEGDAHVDSTTYDGTTPLHIAAGRGSTRLAALLKAAGADPLVENFEPLYDLDDSWENAGEDEGVVPGTTPLDMATSWQVFDILNGKPYEPEFTSDDLLAQGDMKQLAEDVKLQLYKLLEIPDPDKNWATLAQKLGLGILNNAFRLSPAPSKTLMDNYEVSGGTVRELVEALRQMGYTEAIEVIQAASSPVKTTSQAHSLPLSPASTRQQIDELRDSDSVCDSGVETSFRKLSFTESLTSGASLLTLNKMPHDYGQEGPLEGKIX"
## [5,] "MAEDDPYLGRPEQMFHLDPSLTHTIFNPEVFQPQMALPTDGPYLQILEQPKQRGFRFRYVCEGPSHGGLPGASSEKNKKSYPQVKICNYVGPAKVIVQLVTNGKNIHLHAHSLVGKHCEDGICTVTAGPKDMVVGFANLGILHVTKKKVFETLEARMTEACIRGYNPGLLVHPDLAYLQAEGGGDRQLGDREKELIRQAALQQTKEMDLSVVRLMFTAFLPDSTGSFTRRLEPVVSDAIYDSKAPNASNLKIVRMDRTAGCVTGGEEIYLLCDKVQKDDIQIRFYEEEENGGVWEGFGDFSPTDVHRQFAIVFKTPKYKDINITKPASVFVQLRRKSDLETSEPKPFLYYPEIKDKEEVQRKRQKLMPNFSDSFGGGSGAGAGGGGMFGSGGGGGGTGSTGPGYSFPHYGFPTYGGITFHPGTTKSNAGMKHGTMDTESKKDPEGCDKSDDKNTVNLFGKVIETTEQDQEPSEATVGNGEVTLTYATGTKEESAGVQDNLFLEKAMQLAKRHANALFDYAVTGDVKMLLAVQRHLTAVQDENGDSVLHLAIIHLHSQLVRDLLEVTSGLISDDIINMRNDLYQTPLHLAVITKQEDVVEDLLRAGADLSLLDRLGNSVLHLAAKEGHDKVLSILLKHKKAALLLDHPNGDGLNAIHLAMMSNSLPCLLLLVAAGADVNAQEQKSGRTALHLAVEHDNISLAGCLLLEGDAHVDSTTYDGTTPLHIAAGRGSTRLAALLKAAGADPLVENFEPLYDLDDSWENAGEDEGVVPGTTPLDMATSWQVFDILNGKPYEPEFTSDDLLAQGDMKQLAEDVKLQLYKLLEIPDPDKNWATLAQKLGLGILNNAFRLSPAPSKTLMDNYEVSGGTVRELVEALRQMGYTEAIEVIQAASSPVKTTSQAHSLPLSPASTRQQIDELRDSDSVCDSGVETSFRKLSFTESLTSGASLLTLNKMPHDYGQEGPLEGKIXXMAEDDPYLGRPEQMFHLDPSLTHTIFNPEVFQPQMALPTADGPYLQILEQPKQRGFRFRYVCEGPSHGGLPGASSEKNKKSYPQVKICNYVGPAKVIVQLVTNGKNIHLHAHSLVGKHCEDGICTVTAGPKDMVVGFANLGILHVTKKKVFETLEARMTEACIRGYNPGLLVHPDLAYLQAEGGGDRQLGDREKELIRQAALQQTKEMDLSVVRLMFTAFLPDSTGSFTRRLEPVVSDAIYDSKAPNASNLKIVRMDRTAGCVTGGEEIYLLCDKVQKDDIQIRFYEEEENGGVWEGFGDFSPTDVHRQFAIVFKTPKYKDINITKPASVFVQLRRKSDLETSEPKPFLYYPEIKDKEEVQRKRQKLMPNFSDSFGGGSGAGAGGGGMFGSGGGGGGTGSTGPGYSFPHYGFPTYGGITFHPGTTKSNAGMKHGTMDTESKKDPEGCDKSDDKNTVNLFGKVIETTEQDQEPSEATVGNGEVTLTYATGTKEESAGVQDNLFLEKAMQLAKRHANALFDYAVTGDVKMLLAVQRHLTAVQDENGDSVLHLAIIHLHSQLVRDLLEVTSGLISDDIINMRNDLYQTPLHLAVITKQEDVVEDLLRAGADLSLLDRLGNSVLHLAAKEGHDKVLSILLKHKKAALLLDHPNGDGLNAIHLAMMSNSLPCLLLLVAAGADVNAQEQKSGRTALHLAVEHDNISLAGCLLLEGDAHVDSTTYDGTTPLHIAAGRGSTRLAALLKAAGADPLVENFEPLYDLDDSWENAGEDEGVVPGTTPLDMATSWQVFDILNGKPYEPEFTSDDLLAQGDMKQLAEDVKLQLYKLLEIPDPDKNWATLAQKLGLGILNNAFRLSPAPSKTLMDNYEVSGGTVRELVEALRQMGYTEAIEVIQAASSPVKTTSQAHSLPLSPASTRQQIDELRDSDSVCDSGVETSFRKLSFTESLTSGASLLTLNKMPHDYGQEGPLEGKIX"
## [6,] "MAEDDPYLGRPEQMFHLDPSLTHTIFNPEVFQPQMALPTDGPYLQILEQPKQRGFRFRYVCEGPSHGGLPGASSEKNKKSYPQVKICNYVGPAKVIVQLVTNGKNIHLHAHSLVGKHCEDGICTVTAGPKDMVVGFANLGILHVTKKKVFETLEARMTEACIRGYNPGLLVHPDLAYLQAEGGGDRQLGDREKELIRQAALQQTKEMDLSVVRLMFTAFLPDSTGSFTRRLEPVVSDAIYDSKAPNASNLKIVRMDRTAGCVTGGEEIYLLCDKVQKDDIQIRFYEEEENGGVWEGFGDFSPTDVHRQFAIVFKTPKYKDINITKPASVFVQLRRKSDLETSEPKPFLYYPEIKDKEEVQRKRQKLMPNFSDSFGGGSGAGAGGGGMFGSGGGGGGTGSTGPGYSFPHYGFPTYGGITFHPGTTKSNAGMKHGTMDTESKKDPEGCDKSDDKNTVNLFGKVIETTEQDQEPSEATVGNGEVTLTYATGTKEESAGVQDNLFLEKAMQLAKRHANALFDYAVTGDVKMLLAVQRHLTAVQDENGDSVLHLAIIHLHSQLVRDLLEVTSGLISDDIINMRNDLYQTPLHLAVITKQEDVVEDLLRAGADLSLLDRLGNSVLHLAAKEGHDKVLSILLKHKKAALLLDHPNGDGLNAIHLAMMSNSLPCLLLLVAAGADVNAQEQKSGRTALHLAVEHDNISLAGCLLLEGDAHVDSTTYDGTTPLHIAAGRGSTRLAALLKAAGADPLVENFEPLYDLDDSWENAGEDEGVVPGTTPLDMATSWQVFDILNGKPYEPEFTSDDLLAQGDMKQLAEDVKLQLYKLLEIPDPDKNWATLAQKLGLGILNNAFRLSPAPSKTLMDNYEVSGGTVRELVEALRQMGYTEAIEVIQAASSPVKTTSQAHSLPLSPASTRQQIDELRDSDSVCDSGVETSFRKLSFTESLTSGASLLTLNKMPHDYGQEGPLEGKIXXMAEDDPYLGRPEQMFHLDPSLTHTIFNPEVFQPQMALPTADGPYLQILEQPKQRGFRFRYVCEGPSHGGLPGASSEKNKKSYPQVKICNYVGPAKVIVQLVTNGKNIHLHAHSLVGKHCEDGICTVTAGPKDMVVGFANLGILHVTKKKVFETLEARMTEACIRGYNPGLLVHPDLAYLQAEGGGDRQLGDREKELIRQAALQQTKEMDLSVVRLMFTAFLPDSTGSFTRRLEPVVSDAIYDSKAPNASNLKIVRMDRTAGCVTGGEEIYLLCDKVQKDDIQIRFYEEEENGGVWEGFGDFSPTDVHRQFAIVFKTPKYKDINITKPASVFVQLRRKSDLETSEPKPFLYYPEIKDKEEVQRKRQKLMPNFSDSFGGGSGAGAGGGGMFGSGGGGGGTGSTGPGYSFPHYGFPTYGGITFHPGTTKSNAGMKHGTMDTESKKDPEGCDKSDDKNTVNLFGKVIETTEQDQEPSEATVGNGEVTLTYATGTKEESAGVQDNLFLEKAMQLAKRHANALFDYAVTGDVKMLLAVQRHLTAVQDENGDSVLHLAIIHLHSQLVRDLLEVTSGLISDDIINMRNDLYQTPLHLAVITKQEDVVEDLLRAGADLSLLDRLGNSVLHLAAKEGHDKVLSILLKHKKAALLLDHPNGDGLNAIHLAMMSNSLPCLLLLVAAGADVNAQEQKSGRTALHLAVEHDNISLAGCLLLEGDAHVDSTTYDGTTPLHIAAGRGSTRLAALLKAAGADPLVENFEPLYDLDDSWENAGEDEGVVPGTTPLDMATSWQVFDILNGKPYEPEFTSDDLLAQGDMKQLAEDVKLQLYKLLEIPDPDKNWATLAQKLGLGILNNAFRLSPAPSKTLMDNYEVSGGTVRELVEALRQMGYTEAIEVIQAASSPVKTTSQAHSLPLSPASTRQQIDELRDSDSVCDSGVETSFRKLSFTESLTSGASLLTLNKMPHDYGQEGPLEGKIX"
##      Mutant_Peptide                Total_RNA Tumor_RNA_Ratio Tumor_RNA
## [1,] "MAEDDPYLGRPEKMFHLDPSLTHTIFN" "NA"      "NA"            "NA"     
## [2,] "MAEDDPYLGRPEKMFHLDPSLTHTIFN" "NA"      "NA"            "NA"     
## [3,] "MAEDDPYLGRPEKMFHLDPSLTHTIFN" "NA"      "NA"            "NA"     
## [4,] "MAEDDPYLGRPEKMFHLDPSLTHTIFN" "NA"      "NA"            "NA"     
## [5,] "MAEDDPYLGRPEKMFHLDPSLTHTIFN" "NA"      "NA"            "NA"     
## [6,] "MAEDDPYLGRPEKMFHLDPSLTHTIFN" "NA"      "NA"            "NA"     
##      Tumor_RNA_based_on_DNA MutRatio MutRatio_Min MutRatio_Max
## [1,] "NA"                   "NA"     "NA"         "NA"        
## [2,] "NA"                   "NA"     "NA"         "NA"        
## [3,] "NA"                   "NA"     "NA"         "NA"        
## [4,] "NA"                   "NA"     "NA"         "NA"        
## [5,] "NA"                   "NA"     "NA"         "NA"        
## [6,] "NA"                   "NA"     "NA"         "NA"
```

```
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
## 
## [[2]]
##                                                         
## Num_Peptide_Per_NM                                63.000
## Num_Cond_Peptide_Per_NM                           63.000
## Num_Rest_Peptide_Per_NM                            5.000
## Num_Rest_Peptide_Per_NM / Num_Cond_Peptide_Per_NM  0.079
## -logP                                              0.621
## 
## [[3]]
##                                                     MAEDDPYLGRPEKMFHLDPSLTHTIFN-0_atggcagaag-0_1
## Num_Peptide_Per_Pep                                                                       63.000
## Num_Cond_Peptide_Per_Pep                                                                  63.000
## Num_Rest_Peptide_Per_Pep                                                                   5.000
## Num_Rest_Peptide_Per_Pep / Num_Cond_Peptide_Per_Pep                                        0.079
## -logP                                                                                      0.621
```

```
print(head(Result_HLA2_Seq))
```

```
##      HLA                     Pos Gene           Evaluated_Mutant_Peptide_Core
## [1,] "HLA-DPA10103-DPB10201" "5" "0_atggcagaag" "GRPEQMFHL"                  
## [2,] "HLA-DPA10103-DPB10201" "4" "0_atggcagaag" "GRPEQMFHL"                  
## [3,] "HLA-DPA10103-DPB10201" "6" "0_atggcagaag" "EQMFHLILL"                  
## [4,] "DRB1_1302"             "3" "0_atggcagaag" "YLGRPEQMF"                  
## [5,] "DRB1_1302"             "4" "0_atggcagaag" "GRPEQMFHL"                  
## [6,] "DRB1_1302"             "3" "0_atggcagaag" "GRPEQMFHL"
##      Evaluated_Mutant_Peptide Mut_IC50  Mut_Rank Chr NM_ID ReadingFrame
## [1,] "DDPYLGRPEQMFHLI"        "2116.65" "55.00"  "0" ""    "1"         
## [2,] "DPYLGRPEQMFHLIL"        "938.03"  "34.00"  "0" ""    "1"         
## [3,] "PYLGRPEQMFHLILL"        "285.69"  "14.00"  "0" ""    "1"         
## [4,] "DDPYLGRPEQMFHLI"        "5199.08" "90.00"  "0" ""    "1"         
## [5,] "DPYLGRPEQMFHLIL"        "3820.74" "80.00"  "0" ""    "1"         
## [6,] "PYLGRPEQMFHLILL"        "3239.61" "80.00"  "0" ""    "1"         
##      SequenceNumber Chrs                  
## [1,] "1"            "chr4;chr4;chr19;chr4"
## [2,] "1"            "chr4;chr4;chr19;chr4"
## [3,] "1"            "chr4;chr4;chr19;chr4"
## [4,] "1"            "chr4;chr4;chr19;chr4"
## [5,] "1"            "chr4;chr4;chr19;chr4"
## [6,] "1"            "chr4;chr4;chr19;chr4"
##      NM_IDs                                          GeneIDs                 
## [1,] "NM_001165412;NM_003998;NM_005178;NM_001319226" "NFKB1;NFKB1;BCL3;NFKB1"
## [2,] "NM_001165412;NM_003998;NM_005178;NM_001319226" "NFKB1;NFKB1;BCL3;NFKB1"
## [3,] "NM_001165412;NM_003998;NM_005178;NM_001319226" "NFKB1;NFKB1;BCL3;NFKB1"
## [4,] "NM_001165412;NM_003998;NM_005178;NM_001319226" "NFKB1;NFKB1;BCL3;NFKB1"
## [5,] "NM_001165412;NM_003998;NM_005178;NM_001319226" "NFKB1;NFKB1;BCL3;NFKB1"
## [6,] "NM_001165412;NM_003998;NM_005178;NM_001319226" "NFKB1;NFKB1;BCL3;NFKB1"
##      Exon_Starts                             
## [1,] "103422515;103422515;45251964;103423090"
## [2,] "103422515;103422515;45251964;103423090"
## [3,] "103422515;103422515;45251964;103423090"
## [4,] "103422515;103422515;45251964;103423090"
## [5,] "103422515;103422515;45251964;103423090"
## [6,] "103422515;103422515;45251964;103423090"
##      Exon_Ends                                GroupID NumOfPeptides NumOfStops
## [1,] "103538459;103538459;45263301;103538459" "0_1"   "27"          "1"       
## [2,] "103538459;103538459;45263301;103538459" "0_1"   "27"          "1"       
## [3,] "103538459;103538459;45263301;103538459" "0_1"   "27"          "1"       
## [4,] "103538459;103538459;45263301;103538459" "0_1"   "27"          "1"       
## [5,] "103538459;103538459;45263301;103538459" "0_1"   "27"          "1"       
## [6,] "103538459;103538459;45263301;103538459" "0_1"   "27"          "1"       
##      Wt_Peptide                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              
## [1,] "MAEDDPYLGRPEQMFHLDPSLTHTIFNPEVFQPQMALPTDGPYLQILEQPKQRGFRFRYVCEGPSHGGLPGASSEKNKKSYPQVKICNYVGPAKVIVQLVTNGKNIHLHAHSLVGKHCEDGICTVTAGPKDMVVGFANLGILHVTKKKVFETLEARMTEACIRGYNPGLLVHPDLAYLQAEGGGDRQLGDREKELIRQAALQQTKEMDLSVVRLMFTAFLPDSTGSFTRRLEPVVSDAIYDSKAPNASNLKIVRMDRTAGCVTGGEEIYLLCDKVQKDDIQIRFYEEEENGGVWEGFGDFSPTDVHRQFAIVFKTPKYKDINITKPASVFVQLRRKSDLETSEPKPFLYYPEIKDKEEVQRKRQKLMPNFSDSFGGGSGAGAGGGGMFGSGGGGGGTGSTGPGYSFPHYGFPTYGGITFHPGTTKSNAGMKHGTMDTESKKDPEGCDKSDDKNTVNLFGKVIETTEQDQEPSEATVGNGEVTLTYATGTKEESAGVQDNLFLEKAMQLAKRHANALFDYAVTGDVKMLLAVQRHLTAVQDENGDSVLHLAIIHLHSQLVRDLLEVTSGLISDDIINMRNDLYQTPLHLAVITKQEDVVEDLLRAGADLSLLDRLGNSVLHLAAKEGHDKVLSILLKHKKAALLLDHPNGDGLNAIHLAMMSNSLPCLLLLVAAGADVNAQEQKSGRTALHLAVEHDNISLAGCLLLEGDAHVDSTTYDGTTPLHIAAGRGSTRLAALLKAAGADPLVENFEPLYDLDDSWENAGEDEGVVPGTTPLDMATSWQVFDILNGKPYEPEFTSDDLLAQGDMKQLAEDVKLQLYKLLEIPDPDKNWATLAQKLGLGILNNAFRLSPAPSKTLMDNYEVSGGTVRELVEALRQMGYTEAIEVIQAASSPVKTTSQAHSLPLSPASTRQQIDELRDSDSVCDSGVETSFRKLSFTESLTSGASLLTLNKMPHDYGQEGPLEGKIXXMAEDDPYLGRPEQMFHLDPSLTHTIFNPEVFQPQMALPTADGPYLQILEQPKQRGFRFRYVCEGPSHGGLPGASSEKNKKSYPQVKICNYVGPAKVIVQLVTNGKNIHLHAHSLVGKHCEDGICTVTAGPKDMVVGFANLGILHVTKKKVFETLEARMTEACIRGYNPGLLVHPDLAYLQAEGGGDRQLGDREKELIRQAALQQTKEMDLSVVRLMFTAFLPDSTGSFTRRLEPVVSDAIYDSKAPNASNLKIVRMDRTAGCVTGGEEIYLLCDKVQKDDIQIRFYEEEENGGVWEGFGDFSPTDVHRQFAIVFKTPKYKDINITKPASVFVQLRRKSDLETSEPKPFLYYPEIKDKEEVQRKRQKLMPNFSDSFGGGSGAGAGGGGMFGSGGGGGGTGSTGPGYSFPHYGFPTYGGITFHPGTTKSNAGMKHGTMDTESKKDPEGCDKSDDKNTVNLFGKVIETTEQDQEPSEATVGNGEVTLTYATGTKEESAGVQDNLFLEKAMQLAKRHANALFDYAVTGDVKMLLAVQRHLTAVQDENGDSVLHLAIIHLHSQLVRDLLEVTSGLISDDIINMRNDLYQTPLHLAVITKQEDVVEDLLRAGADLSLLDRLGNSVLHLAAKEGHDKVLSILLKHKKAALLLDHPNGDGLNAIHLAMMSNSLPCLLLLVAAGADVNAQEQKSGRTALHLAVEHDNISLAGCLLLEGDAHVDSTTYDGTTPLHIAAGRGSTRLAALLKAAGADPLVENFEPLYDLDDSWENAGEDEGVVPGTTPLDMATSWQVFDILNGKPYEPEFTSDDLLAQGDMKQLAEDVKLQLYKLLEIPDPDKNWATLAQKLGLGILNNAFRLSPAPSKTLMDNYEVSGGTVRELVEALRQMGYTEAIEVIQAASSPVKTTSQAHSLPLSPASTRQQIDELRDSDSVCDSGVETSFRKLSFTESLTSGASLLTLNKMPHDYGQEGPLEGKIXXMPRCPAGAMDEGPVDLRTRPKAAGLPGAALPLRKRPLRAPSPEPAAPRGAAGLVVPLDPLRGGCDLPAVPGPPHGLARPEALYYPGALLPLYPTRAMGSPFPLVNLPTPLYPMMCPMEHPLSADIAMATRADEDGDTPLHIAVVQGNLPAVHRLVNLFQQGGRELDIYNNLRQTPLHLAVITTLPSVVRLLVTAGASPMALDRHGQTAAHLACEHRSPTCLRALLDSAAPGTLDLEARNYDGLTALHVAVNTECQETVQLLLERGADIDAVDIKSGRSPLIHAVENNSLSMVQLLLQHGANVNAQMYSGSSALHSASGRGLLPLVRTLVRSGADSSLKNCHNDTPLMVARSRRVIDILRGKATRPASTSQPDPSPDRSANTSPESSSRLSSNGLLSASPSSSPSQSPPRDPPGFPMAPPNFFLPSPSPPAFLPFAGVLRGPGRPVPPSPAPGGSXXMAEDDPYLGRPEQMFHLDPSLTHTIFNPEVFQPQMALPTDGPYLQILEQPKQRGFRFRYVCEGPSHGGLPGASSEKNKKSYPQVKICNYVGPAKVIVQLVTNGKNIHLHAHSLVGKHCEDGICTVTAGPKDMVVGFANLGILHVTKKKVFETLEARMTEACIRGYNPGLLVHPDLAYLQAEGGGDRQLGDREKELIRQAALQQTKEMDLSVVRLMFTAFLPDSTGSFTRRLEPVVSDAIYDSKAPNASNLKIVRMDRTAGCVTGGEEIYLLCDKVQKDDIQIRFYEEEENGGVWEGFGDFSPTDVHRQFAIVFKTPKYKDINITKPASVFVQLRRKSDLETSEPKPFLYYPEIKDKEEVQRKRQKLMPNFSDSFGGGSGAGAGGGGMFGSGGGGGGTGSTGPGYSFPHYGFPTYGGITFHPGTTKSNAGMKHGTMDTESKKDPEGCDKSDDKNTVNLFGKVIETTEQDQEPSEATVGNGEVTLTYATGTKEESAGVQDNLFLEKAMQLAKRHANALFDYAVTGDVKMLLAVQRHLTAVQDENGDSVLHLAIIHLHSQLVRDLLEVTSGLISDDIINMRNDLYQTPLHLAVITKQEDVVEDLLRAGADLSLLDRLGNSVLHLAAKEGHDKVLSILLKHKKAALLLDHPNGDGLNAIHLAMMSNSLPCLLLLVAAGADVNAQEQKSGRTALHLAVEHDNISLAGCLLLEGDAHVDSTTYDGTTPLHIAAGRGSTRLAALLKAAGADPLVENFEPLYDLDDSWENAGEDEGVVPGTTPLDMATSWQVFDILNGKPYEPEFTSDDLLAQGDMKQLAEDVKLQLYKLLEIPDPDKNWATLAQKLGLGILNNAFRLSPAPSKTLMDNYEVSGGTVRELVEALRQMGYTEAIEVIQAASSPVKTTSQAHSLPLSPASTRQQIDELRDSDSVCDSGVETSFRKLSFTESLTSGASLLTLNKMPHDYGQEGPLEGKIX"
## [2,] "MAEDDPYLGRPEQMFHLDPSLTHTIFNPEVFQPQMALPTDGPYLQILEQPKQRGFRFRYVCEGPSHGGLPGASSEKNKKSYPQVKICNYVGPAKVIVQLVTNGKNIHLHAHSLVGKHCEDGICTVTAGPKDMVVGFANLGILHVTKKKVFETLEARMTEACIRGYNPGLLVHPDLAYLQAEGGGDRQLGDREKELIRQAALQQTKEMDLSVVRLMFTAFLPDSTGSFTRRLEPVVSDAIYDSKAPNASNLKIVRMDRTAGCVTGGEEIYLLCDKVQKDDIQIRFYEEEENGGVWEGFGDFSPTDVHRQFAIVFKTPKYKDINITKPASVFVQLRRKSDLETSEPKPFLYYPEIKDKEEVQRKRQKLMPNFSDSFGGGSGAGAGGGGMFGSGGGGGGTGSTGPGYSFPHYGFPTYGGITFHPGTTKSNAGMKHGTMDTESKKDPEGCDKSDDKNTVNLFGKVIETTEQDQEPSEATVGNGEVTLTYATGTKEESAGVQDNLFLEKAMQLAKRHANALFDYAVTGDVKMLLAVQRHLTAVQDENGDSVLHLAIIHLHSQLVRDLLEVTSGLISDDIINMRNDLYQTPLHLAVITKQEDVVEDLLRAGADLSLLDRLGNSVLHLAAKEGHDKVLSILLKHKKAALLLDHPNGDGLNAIHLAMMSNSLPCLLLLVAAGADVNAQEQKSGRTALHLAVEHDNISLAGCLLLEGDAHVDSTTYDGTTPLHIAAGRGSTRLAALLKAAGADPLVENFEPLYDLDDSWENAGEDEGVVPGTTPLDMATSWQVFDILNGKPYEPEFTSDDLLAQGDMKQLAEDVKLQLYKLLEIPDPDKNWATLAQKLGLGILNNAFRLSPAPSKTLMDNYEVSGGTVRELVEALRQMGYTEAIEVIQAASSPVKTTSQAHSLPLSPASTRQQIDELRDSDSVCDSGVETSFRKLSFTESLTSGASLLTLNKMPHDYGQEGPLEGKIXXMAEDDPYLGRPEQMFHLDPSLTHTIFNPEVFQPQMALPTADGPYLQILEQPKQRGFRFRYVCEGPSHGGLPGASSEKNKKSYPQVKICNYVGPAKVIVQLVTNGKNIHLHAHSLVGKHCEDGICTVTAGPKDMVVGFANLGILHVTKKKVFETLEARMTEACIRGYNPGLLVHPDLAYLQAEGGGDRQLGDREKELIRQAALQQTKEMDLSVVRLMFTAFLPDSTGSFTRRLEPVVSDAIYDSKAPNASNLKIVRMDRTAGCVTGGEEIYLLCDKVQKDDIQIRFYEEEENGGVWEGFGDFSPTDVHRQFAIVFKTPKYKDINITKPASVFVQLRRKSDLETSEPKPFLYYPEIKDKEEVQRKRQKLMPNFSDSFGGGSGAGAGGGGMFGSGGGGGGTGSTGPGYSFPHYGFPTYGGITFHPGTTKSNAGMKHGTMDTESKKDPEGCDKSDDKNTVNLFGKVIETTEQDQEPSEATVGNGEVTLTYATGTKEESAGVQDNLFLEKAMQLAKRHANALFDYAVTGDVKMLLAVQRHLTAVQDENGDSVLHLAIIHLHSQLVRDLLEVTSGLISDDIINMRNDLYQTPLHLAVITKQEDVVEDLLRAGADLSLLDRLGNSVLHLAAKEGHDKVLSILLKHKKAALLLDHPNGDGLNAIHLAMMSNSLPCLLLLVAAGADVNAQEQKSGRTALHLAVEHDNISLAGCLLLEGDAHVDSTTYDGTTPLHIAAGRGSTRLAALLKAAGADPLVENFEPLYDLDDSWENAGEDEGVVPGTTPLDMATSWQVFDILNGKPYEPEFTSDDLLAQGDMKQLAEDVKLQLYKLLEIPDPDKNWATLAQKLGLGILNNAFRLSPAPSKTLMDNYEVSGGTVRELVEALRQMGYTEAIEVIQAASSPVKTTSQAHSLPLSPASTRQQIDELRDSDSVCDSGVETSFRKLSFTESLTSGASLLTLNKMPHDYGQEGPLEGKIXXMPRCPAGAMDEGPVDLRTRPKAAGLPGAALPLRKRPLRAPSPEPAAPRGAAGLVVPLDPLRGGCDLPAVPGPPHGLARPEALYYPGALLPLYPTRAMGSPFPLVNLPTPLYPMMCPMEHPLSADIAMATRADEDGDTPLHIAVVQGNLPAVHRLVNLFQQGGRELDIYNNLRQTPLHLAVITTLPSVVRLLVTAGASPMALDRHGQTAAHLACEHRSPTCLRALLDSAAPGTLDLEARNYDGLTALHVAVNTECQETVQLLLERGADIDAVDIKSGRSPLIHAVENNSLSMVQLLLQHGANVNAQMYSGSSALHSASGRGLLPLVRTLVRSGADSSLKNCHNDTPLMVARSRRVIDILRGKATRPASTSQPDPSPDRSANTSPESSSRLSSNGLLSASPSSSPSQSPPRDPPGFPMAPPNFFLPSPSPPAFLPFAGVLRGPGRPVPPSPAPGGSXXMAEDDPYLGRPEQMFHLDPSLTHTIFNPEVFQPQMALPTDGPYLQILEQPKQRGFRFRYVCEGPSHGGLPGASSEKNKKSYPQVKICNYVGPAKVIVQLVTNGKNIHLHAHSLVGKHCEDGICTVTAGPKDMVVGFANLGILHVTKKKVFETLEARMTEACIRGYNPGLLVHPDLAYLQAEGGGDRQLGDREKELIRQAALQQTKEMDLSVVRLMFTAFLPDSTGSFTRRLEPVVSDAIYDSKAPNASNLKIVRMDRTAGCVTGGEEIYLLCDKVQKDDIQIRFYEEEENGGVWEGFGDFSPTDVHRQFAIVFKTPKYKDINITKPASVFVQLRRKSDLETSEPKPFLYYPEIKDKEEVQRKRQKLMPNFSDSFGGGSGAGAGGGGMFGSGGGGGGTGSTGPGYSFPHYGFPTYGGITFHPGTTKSNAGMKHGTMDTESKKDPEGCDKSDDKNTVNLFGKVIETTEQDQEPSEATVGNGEVTLTYATGTKEESAGVQDNLFLEKAMQLAKRHANALFDYAVTGDVKMLLAVQRHLTAVQDENGDSVLHLAIIHLHSQLVRDLLEVTSGLISDDIINMRNDLYQTPLHLAVITKQEDVVEDLLRAGADLSLLDRLGNSVLHLAAKEGHDKVLSILLKHKKAALLLDHPNGDGLNAIHLAMMSNSLPCLLLLVAAGADVNAQEQKSGRTALHLAVEHDNISLAGCLLLEGDAHVDSTTYDGTTPLHIAAGRGSTRLAALLKAAGADPLVENFEPLYDLDDSWENAGEDEGVVPGTTPLDMATSWQVFDILNGKPYEPEFTSDDLLAQGDMKQLAEDVKLQLYKLLEIPDPDKNWATLAQKLGLGILNNAFRLSPAPSKTLMDNYEVSGGTVRELVEALRQMGYTEAIEVIQAASSPVKTTSQAHSLPLSPASTRQQIDELRDSDSVCDSGVETSFRKLSFTESLTSGASLLTLNKMPHDYGQEGPLEGKIX"
## [3,] "MAEDDPYLGRPEQMFHLDPSLTHTIFNPEVFQPQMALPTDGPYLQILEQPKQRGFRFRYVCEGPSHGGLPGASSEKNKKSYPQVKICNYVGPAKVIVQLVTNGKNIHLHAHSLVGKHCEDGICTVTAGPKDMVVGFANLGILHVTKKKVFETLEARMTEACIRGYNPGLLVHPDLAYLQAEGGGDRQLGDREKELIRQAALQQTKEMDLSVVRLMFTAFLPDSTGSFTRRLEPVVSDAIYDSKAPNASNLKIVRMDRTAGCVTGGEEIYLLCDKVQKDDIQIRFYEEEENGGVWEGFGDFSPTDVHRQFAIVFKTPKYKDINITKPASVFVQLRRKSDLETSEPKPFLYYPEIKDKEEVQRKRQKLMPNFSDSFGGGSGAGAGGGGMFGSGGGGGGTGSTGPGYSFPHYGFPTYGGITFHPGTTKSNAGMKHGTMDTESKKDPEGCDKSDDKNTVNLFGKVIETTEQDQEPSEATVGNGEVTLTYATGTKEESAGVQDNLFLEKAMQLAKRHANALFDYAVTGDVKMLLAVQRHLTAVQDENGDSVLHLAIIHLHSQLVRDLLEVTSGLISDDIINMRNDLYQTPLHLAVITKQEDVVEDLLRAGADLSLLDRLGNSVLHLAAKEGHDKVLSILLKHKKAALLLDHPNGDGLNAIHLAMMSNSLPCLLLLVAAGADVNAQEQKSGRTALHLAVEHDNISLAGCLLLEGDAHVDSTTYDGTTPLHIAAGRGSTRLAALLKAAGADPLVENFEPLYDLDDSWENAGEDEGVVPGTTPLDMATSWQVFDILNGKPYEPEFTSDDLLAQGDMKQLAEDVKLQLYKLLEIPDPDKNWATLAQKLGLGILNNAFRLSPAPSKTLMDNYEVSGGTVRELVEALRQMGYTEAIEVIQAASSPVKTTSQAHSLPLSPASTRQQIDELRDSDSVCDSGVETSFRKLSFTESLTSGASLLTLNKMPHDYGQEGPLEGKIXXMAEDDPYLGRPEQMFHLDPSLTHTIFNPEVFQPQMALPTADGPYLQILEQPKQRGFRFRYVCEGPSHGGLPGASSEKNKKSYPQVKICNYVGPAKVIVQLVTNGKNIHLHAHSLVGKHCEDGICTVTAGPKDMVVGFANLGILHVTKKKVFETLEARMTEACIRGYNPGLLVHPDLAYLQAEGGGDRQLGDREKELIRQAALQQTKEMDLSVVRLMFTAFLPDSTGSFTRRLEPVVSDAIYDSKAPNASNLKIVRMDRTAGCVTGGEEIYLLCDKVQKDDIQIRFYEEEENGGVWEGFGDFSPTDVHRQFAIVFKTPKYKDINITKPASVFVQLRRKSDLETSEPKPFLYYPEIKDKEEVQRKRQKLMPNFSDSFGGGSGAGAGGGGMFGSGGGGGGTGSTGPGYSFPHYGFPTYGGITFHPGTTKSNAGMKHGTMDTESKKDPEGCDKSDDKNTVNLFGKVIETTEQDQEPSEATVGNGEVTLTYATGTKEESAGVQDNLFLEKAMQLAKRHANALFDYAVTGDVKMLLAVQRHLTAVQDENGDSVLHLAIIHLHSQLVRDLLEVTSGLISDDIINMRNDLYQTPLHLAVITKQEDVVEDLLRAGADLSLLDRLGNSVLHLAAKEGHDKVLSILLKHKKAALLLDHPNGDGLNAIHLAMMSNSLPCLLLLVAAGADVNAQEQKSGRTALHLAVEHDNISLAGCLLLEGDAHVDSTTYDGTTPLHIAAGRGSTRLAALLKAAGADPLVENFEPLYDLDDSWENAGEDEGVVPGTTPLDMATSWQVFDILNGKPYEPEFTSDDLLAQGDMKQLAEDVKLQLYKLLEIPDPDKNWATLAQKLGLGILNNAFRLSPAPSKTLMDNYEVSGGTVRELVEALRQMGYTEAIEVIQAASSPVKTTSQAHSLPLSPASTRQQIDELRDSDSVCDSGVETSFRKLSFTESLTSGASLLTLNKMPHDYGQEGPLEGKIXXMPRCPAGAMDEGPVDLRTRPKAAGLPGAALPLRKRPLRAPSPEPAAPRGAAGLVVPLDPLRGGCDLPAVPGPPHGLARPEALYYPGALLPLYPTRAMGSPFPLVNLPTPLYPMMCPMEHPLSADIAMATRADEDGDTPLHIAVVQGNLPAVHRLVNLFQQGGRELDIYNNLRQTPLHLAVITTLPSVVRLLVTAGASPMALDRHGQTAAHLACEHRSPTCLRALLDSAAPGTLDLEARNYDGLTALHVAVNTECQETVQLLLERGADIDAVDIKSGRSPLIHAVENNSLSMVQLLLQHGANVNAQMYSGSSALHSASGRGLLPLVRTLVRSGADSSLKNCHNDTPLMVARSRRVIDILRGKATRPASTSQPDPSPDRSANTSPESSSRLSSNGLLSASPSSSPSQSPPRDPPGFPMAPPNFFLPSPSPPAFLPFAGVLRGPGRPVPPSPAPGGSXXMAEDDPYLGRPEQMFHLDPSLTHTIFNPEVFQPQMALPTDGPYLQILEQPKQRGFRFRYVCEGPSHGGLPGASSEKNKKSYPQVKICNYVGPAKVIVQLVTNGKNIHLHAHSLVGKHCEDGICTVTAGPKDMVVGFANLGILHVTKKKVFETLEARMTEACIRGYNPGLLVHPDLAYLQAEGGGDRQLGDREKELIRQAALQQTKEMDLSVVRLMFTAFLPDSTGSFTRRLEPVVSDAIYDSKAPNASNLKIVRMDRTAGCVTGGEEIYLLCDKVQKDDIQIRFYEEEENGGVWEGFGDFSPTDVHRQFAIVFKTPKYKDINITKPASVFVQLRRKSDLETSEPKPFLYYPEIKDKEEVQRKRQKLMPNFSDSFGGGSGAGAGGGGMFGSGGGGGGTGSTGPGYSFPHYGFPTYGGITFHPGTTKSNAGMKHGTMDTESKKDPEGCDKSDDKNTVNLFGKVIETTEQDQEPSEATVGNGEVTLTYATGTKEESAGVQDNLFLEKAMQLAKRHANALFDYAVTGDVKMLLAVQRHLTAVQDENGDSVLHLAIIHLHSQLVRDLLEVTSGLISDDIINMRNDLYQTPLHLAVITKQEDVVEDLLRAGADLSLLDRLGNSVLHLAAKEGHDKVLSILLKHKKAALLLDHPNGDGLNAIHLAMMSNSLPCLLLLVAAGADVNAQEQKSGRTALHLAVEHDNISLAGCLLLEGDAHVDSTTYDGTTPLHIAAGRGSTRLAALLKAAGADPLVENFEPLYDLDDSWENAGEDEGVVPGTTPLDMATSWQVFDILNGKPYEPEFTSDDLLAQGDMKQLAEDVKLQLYKLLEIPDPDKNWATLAQKLGLGILNNAFRLSPAPSKTLMDNYEVSGGTVRELVEALRQMGYTEAIEVIQAASSPVKTTSQAHSLPLSPASTRQQIDELRDSDSVCDSGVETSFRKLSFTESLTSGASLLTLNKMPHDYGQEGPLEGKIX"
## [4,] "MAEDDPYLGRPEQMFHLDPSLTHTIFNPEVFQPQMALPTDGPYLQILEQPKQRGFRFRYVCEGPSHGGLPGASSEKNKKSYPQVKICNYVGPAKVIVQLVTNGKNIHLHAHSLVGKHCEDGICTVTAGPKDMVVGFANLGILHVTKKKVFETLEARMTEACIRGYNPGLLVHPDLAYLQAEGGGDRQLGDREKELIRQAALQQTKEMDLSVVRLMFTAFLPDSTGSFTRRLEPVVSDAIYDSKAPNASNLKIVRMDRTAGCVTGGEEIYLLCDKVQKDDIQIRFYEEEENGGVWEGFGDFSPTDVHRQFAIVFKTPKYKDINITKPASVFVQLRRKSDLETSEPKPFLYYPEIKDKEEVQRKRQKLMPNFSDSFGGGSGAGAGGGGMFGSGGGGGGTGSTGPGYSFPHYGFPTYGGITFHPGTTKSNAGMKHGTMDTESKKDPEGCDKSDDKNTVNLFGKVIETTEQDQEPSEATVGNGEVTLTYATGTKEESAGVQDNLFLEKAMQLAKRHANALFDYAVTGDVKMLLAVQRHLTAVQDENGDSVLHLAIIHLHSQLVRDLLEVTSGLISDDIINMRNDLYQTPLHLAVITKQEDVVEDLLRAGADLSLLDRLGNSVLHLAAKEGHDKVLSILLKHKKAALLLDHPNGDGLNAIHLAMMSNSLPCLLLLVAAGADVNAQEQKSGRTALHLAVEHDNISLAGCLLLEGDAHVDSTTYDGTTPLHIAAGRGSTRLAALLKAAGADPLVENFEPLYDLDDSWENAGEDEGVVPGTTPLDMATSWQVFDILNGKPYEPEFTSDDLLAQGDMKQLAEDVKLQLYKLLEIPDPDKNWATLAQKLGLGILNNAFRLSPAPSKTLMDNYEVSGGTVRELVEALRQMGYTEAIEVIQAASSPVKTTSQAHSLPLSPASTRQQIDELRDSDSVCDSGVETSFRKLSFTESLTSGASLLTLNKMPHDYGQEGPLEGKIXXMAEDDPYLGRPEQMFHLDPSLTHTIFNPEVFQPQMALPTADGPYLQILEQPKQRGFRFRYVCEGPSHGGLPGASSEKNKKSYPQVKICNYVGPAKVIVQLVTNGKNIHLHAHSLVGKHCEDGICTVTAGPKDMVVGFANLGILHVTKKKVFETLEARMTEACIRGYNPGLLVHPDLAYLQAEGGGDRQLGDREKELIRQAALQQTKEMDLSVVRLMFTAFLPDSTGSFTRRLEPVVSDAIYDSKAPNASNLKIVRMDRTAGCVTGGEEIYLLCDKVQKDDIQIRFYEEEENGGVWEGFGDFSPTDVHRQFAIVFKTPKYKDINITKPASVFVQLRRKSDLETSEPKPFLYYPEIKDKEEVQRKRQKLMPNFSDSFGGGSGAGAGGGGMFGSGGGGGGTGSTGPGYSFPHYGFPTYGGITFHPGTTKSNAGMKHGTMDTESKKDPEGCDKSDDKNTVNLFGKVIETTEQDQEPSEATVGNGEVTLTYATGTKEESAGVQDNLFLEKAMQLAKRHANALFDYAVTGDVKMLLAVQRHLTAVQDENGDSVLHLAIIHLHSQLVRDLLEVTSGLISDDIINMRNDLYQTPLHLAVITKQEDVVEDLLRAGADLSLLDRLGNSVLHLAAKEGHDKVLSILLKHKKAALLLDHPNGDGLNAIHLAMMSNSLPCLLLLVAAGADVNAQEQKSGRTALHLAVEHDNISLAGCLLLEGDAHVDSTTYDGTTPLHIAAGRGSTRLAALLKAAGADPLVENFEPLYDLDDSWENAGEDEGVVPGTTPLDMATSWQVFDILNGKPYEPEFTSDDLLAQGDMKQLAEDVKLQLYKLLEIPDPDKNWATLAQKLGLGILNNAFRLSPAPSKTLMDNYEVSGGTVRELVEALRQMGYTEAIEVIQAASSPVKTTSQAHSLPLSPASTRQQIDELRDSDSVCDSGVETSFRKLSFTESLTSGASLLTLNKMPHDYGQEGPLEGKIXXMPRCPAGAMDEGPVDLRTRPKAAGLPGAALPLRKRPLRAPSPEPAAPRGAAGLVVPLDPLRGGCDLPAVPGPPHGLARPEALYYPGALLPLYPTRAMGSPFPLVNLPTPLYPMMCPMEHPLSADIAMATRADEDGDTPLHIAVVQGNLPAVHRLVNLFQQGGRELDIYNNLRQTPLHLAVITTLPSVVRLLVTAGASPMALDRHGQTAAHLACEHRSPTCLRALLDSAAPGTLDLEARNYDGLTALHVAVNTECQETVQLLLERGADIDAVDIKSGRSPLIHAVENNSLSMVQLLLQHGANVNAQMYSGSSALHSASGRGLLPLVRTLVRSGADSSLKNCHNDTPLMVARSRRVIDILRGKATRPASTSQPDPSPDRSANTSPESSSRLSSNGLLSASPSSSPSQSPPRDPPGFPMAPPNFFLPSPSPPAFLPFAGVLRGPGRPVPPSPAPGGSXXMAEDDPYLGRPEQMFHLDPSLTHTIFNPEVFQPQMALPTDGPYLQILEQPKQRGFRFRYVCEGPSHGGLPGASSEKNKKSYPQVKICNYVGPAKVIVQLVTNGKNIHLHAHSLVGKHCEDGICTVTAGPKDMVVGFANLGILHVTKKKVFETLEARMTEACIRGYNPGLLVHPDLAYLQAEGGGDRQLGDREKELIRQAALQQTKEMDLSVVRLMFTAFLPDSTGSFTRRLEPVVSDAIYDSKAPNASNLKIVRMDRTAGCVTGGEEIYLLCDKVQKDDIQIRFYEEEENGGVWEGFGDFSPTDVHRQFAIVFKTPKYKDINITKPASVFVQLRRKSDLETSEPKPFLYYPEIKDKEEVQRKRQKLMPNFSDSFGGGSGAGAGGGGMFGSGGGGGGTGSTGPGYSFPHYGFPTYGGITFHPGTTKSNAGMKHGTMDTESKKDPEGCDKSDDKNTVNLFGKVIETTEQDQEPSEATVGNGEVTLTYATGTKEESAGVQDNLFLEKAMQLAKRHANALFDYAVTGDVKMLLAVQRHLTAVQDENGDSVLHLAIIHLHSQLVRDLLEVTSGLISDDIINMRNDLYQTPLHLAVITKQEDVVEDLLRAGADLSLLDRLGNSVLHLAAKEGHDKVLSILLKHKKAALLLDHPNGDGLNAIHLAMMSNSLPCLLLLVAAGADVNAQEQKSGRTALHLAVEHDNISLAGCLLLEGDAHVDSTTYDGTTPLHIAAGRGSTRLAALLKAAGADPLVENFEPLYDLDDSWENAGEDEGVVPGTTPLDMATSWQVFDILNGKPYEPEFTSDDLLAQGDMKQLAEDVKLQLYKLLEIPDPDKNWATLAQKLGLGILNNAFRLSPAPSKTLMDNYEVSGGTVRELVEALRQMGYTEAIEVIQAASSPVKTTSQAHSLPLSPASTRQQIDELRDSDSVCDSGVETSFRKLSFTESLTSGASLLTLNKMPHDYGQEGPLEGKIX"
## [5,] "MAEDDPYLGRPEQMFHLDPSLTHTIFNPEVFQPQMALPTDGPYLQILEQPKQRGFRFRYVCEGPSHGGLPGASSEKNKKSYPQVKICNYVGPAKVIVQLVTNGKNIHLHAHSLVGKHCEDGICTVTAGPKDMVVGFANLGILHVTKKKVFETLEARMTEACIRGYNPGLLVHPDLAYLQAEGGGDRQLGDREKELIRQAALQQTKEMDLSVVRLMFTAFLPDSTGSFTRRLEPVVSDAIYDSKAPNASNLKIVRMDRTAGCVTGGEEIYLLCDKVQKDDIQIRFYEEEENGGVWEGFGDFSPTDVHRQFAIVFKTPKYKDINITKPASVFVQLRRKSDLETSEPKPFLYYPEIKDKEEVQRKRQKLMPNFSDSFGGGSGAGAGGGGMFGSGGGGGGTGSTGPGYSFPHYGFPTYGGITFHPGTTKSNAGMKHGTMDTESKKDPEGCDKSDDKNTVNLFGKVIETTEQDQEPSEATVGNGEVTLTYATGTKEESAGVQDNLFLEKAMQLAKRHANALFDYAVTGDVKMLLAVQRHLTAVQDENGDSVLHLAIIHLHSQLVRDLLEVTSGLISDDIINMRNDLYQTPLHLAVITKQEDVVEDLLRAGADLSLLDRLGNSVLHLAAKEGHDKVLSILLKHKKAALLLDHPNGDGLNAIHLAMMSNSLPCLLLLVAAGADVNAQEQKSGRTALHLAVEHDNISLAGCLLLEGDAHVDSTTYDGTTPLHIAAGRGSTRLAALLKAAGADPLVENFEPLYDLDDSWENAGEDEGVVPGTTPLDMATSWQVFDILNGKPYEPEFTSDDLLAQGDMKQLAEDVKLQLYKLLEIPDPDKNWATLAQKLGLGILNNAFRLSPAPSKTLMDNYEVSGGTVRELVEALRQMGYTEAIEVIQAASSPVKTTSQAHSLPLSPASTRQQIDELRDSDSVCDSGVETSFRKLSFTESLTSGASLLTLNKMPHDYGQEGPLEGKIXXMAEDDPYLGRPEQMFHLDPSLTHTIFNPEVFQPQMALPTADGPYLQILEQPKQRGFRFRYVCEGPSHGGLPGASSEKNKKSYPQVKICNYVGPAKVIVQLVTNGKNIHLHAHSLVGKHCEDGICTVTAGPKDMVVGFANLGILHVTKKKVFETLEARMTEACIRGYNPGLLVHPDLAYLQAEGGGDRQLGDREKELIRQAALQQTKEMDLSVVRLMFTAFLPDSTGSFTRRLEPVVSDAIYDSKAPNASNLKIVRMDRTAGCVTGGEEIYLLCDKVQKDDIQIRFYEEEENGGVWEGFGDFSPTDVHRQFAIVFKTPKYKDINITKPASVFVQLRRKSDLETSEPKPFLYYPEIKDKEEVQRKRQKLMPNFSDSFGGGSGAGAGGGGMFGSGGGGGGTGSTGPGYSFPHYGFPTYGGITFHPGTTKSNAGMKHGTMDTESKKDPEGCDKSDDKNTVNLFGKVIETTEQDQEPSEATVGNGEVTLTYATGTKEESAGVQDNLFLEKAMQLAKRHANALFDYAVTGDVKMLLAVQRHLTAVQDENGDSVLHLAIIHLHSQLVRDLLEVTSGLISDDIINMRNDLYQTPLHLAVITKQEDVVEDLLRAGADLSLLDRLGNSVLHLAAKEGHDKVLSILLKHKKAALLLDHPNGDGLNAIHLAMMSNSLPCLLLLVAAGADVNAQEQKSGRTALHLAVEHDNISLAGCLLLEGDAHVDSTTYDGTTPLHIAAGRGSTRLAALLKAAGADPLVENFEPLYDLDDSWENAGEDEGVVPGTTPLDMATSWQVFDILNGKPYEPEFTSDDLLAQGDMKQLAEDVKLQLYKLLEIPDPDKNWATLAQKLGLGILNNAFRLSPAPSKTLMDNYEVSGGTVRELVEALRQMGYTEAIEVIQAASSPVKTTSQAHSLPLSPASTRQQIDELRDSDSVCDSGVETSFRKLSFTESLTSGASLLTLNKMPHDYGQEGPLEGKIXXMPRCPAGAMDEGPVDLRTRPKAAGLPGAALPLRKRPLRAPSPEPAAPRGAAGLVVPLDPLRGGCDLPAVPGPPHGLARPEALYYPGALLPLYPTRAMGSPFPLVNLPTPLYPMMCPMEHPLSADIAMATRADEDGDTPLHIAVVQGNLPAVHRLVNLFQQGGRELDIYNNLRQTPLHLAVITTLPSVVRLLVTAGASPMALDRHGQTAAHLACEHRSPTCLRALLDSAAPGTLDLEARNYDGLTALHVAVNTECQETVQLLLERGADIDAVDIKSGRSPLIHAVENNSLSMVQLLLQHGANVNAQMYSGSSALHSASGRGLLPLVRTLVRSGADSSLKNCHNDTPLMVARSRRVIDILRGKATRPASTSQPDPSPDRSANTSPESSSRLSSNGLLSASPSSSPSQSPPRDPPGFPMAPPNFFLPSPSPPAFLPFAGVLRGPGRPVPPSPAPGGSXXMAEDDPYLGRPEQMFHLDPSLTHTIFNPEVFQPQMALPTDGPYLQILEQPKQRGFRFRYVCEGPSHGGLPGASSEKNKKSYPQVKICNYVGPAKVIVQLVTNGKNIHLHAHSLVGKHCEDGICTVTAGPKDMVVGFANLGILHVTKKKVFETLEARMTEACIRGYNPGLLVHPDLAYLQAEGGGDRQLGDREKELIRQAALQQTKEMDLSVVRLMFTAFLPDSTGSFTRRLEPVVSDAIYDSKAPNASNLKIVRMDRTAGCVTGGEEIYLLCDKVQKDDIQIRFYEEEENGGVWEGFGDFSPTDVHRQFAIVFKTPKYKDINITKPASVFVQLRRKSDLETSEPKPFLYYPEIKDKEEVQRKRQKLMPNFSDSFGGGSGAGAGGGGMFGSGGGGGGTGSTGPGYSFPHYGFPTYGGITFHPGTTKSNAGMKHGTMDTESKKDPEGCDKSDDKNTVNLFGKVIETTEQDQEPSEATVGNGEVTLTYATGTKEESAGVQDNLFLEKAMQLAKRHANALFDYAVTGDVKMLLAVQRHLTAVQDENGDSVLHLAIIHLHSQLVRDLLEVTSGLISDDIINMRNDLYQTPLHLAVITKQEDVVEDLLRAGADLSLLDRLGNSVLHLAAKEGHDKVLSILLKHKKAALLLDHPNGDGLNAIHLAMMSNSLPCLLLLVAAGADVNAQEQKSGRTALHLAVEHDNISLAGCLLLEGDAHVDSTTYDGTTPLHIAAGRGSTRLAALLKAAGADPLVENFEPLYDLDDSWENAGEDEGVVPGTTPLDMATSWQVFDILNGKPYEPEFTSDDLLAQGDMKQLAEDVKLQLYKLLEIPDPDKNWATLAQKLGLGILNNAFRLSPAPSKTLMDNYEVSGGTVRELVEALRQMGYTEAIEVIQAASSPVKTTSQAHSLPLSPASTRQQIDELRDSDSVCDSGVETSFRKLSFTESLTSGASLLTLNKMPHDYGQEGPLEGKIX"
## [6,] "MAEDDPYLGRPEQMFHLDPSLTHTIFNPEVFQPQMALPTDGPYLQILEQPKQRGFRFRYVCEGPSHGGLPGASSEKNKKSYPQVKICNYVGPAKVIVQLVTNGKNIHLHAHSLVGKHCEDGICTVTAGPKDMVVGFANLGILHVTKKKVFETLEARMTEACIRGYNPGLLVHPDLAYLQAEGGGDRQLGDREKELIRQAALQQTKEMDLSVVRLMFTAFLPDSTGSFTRRLEPVVSDAIYDSKAPNASNLKIVRMDRTAGCVTGGEEIYLLCDKVQKDDIQIRFYEEEENGGVWEGFGDFSPTDVHRQFAIVFKTPKYKDINITKPASVFVQLRRKSDLETSEPKPFLYYPEIKDKEEVQRKRQKLMPNFSDSFGGGSGAGAGGGGMFGSGGGGGGTGSTGPGYSFPHYGFPTYGGITFHPGTTKSNAGMKHGTMDTESKKDPEGCDKSDDKNTVNLFGKVIETTEQDQEPSEATVGNGEVTLTYATGTKEESAGVQDNLFLEKAMQLAKRHANALFDYAVTGDVKMLLAVQRHLTAVQDENGDSVLHLAIIHLHSQLVRDLLEVTSGLISDDIINMRNDLYQTPLHLAVITKQEDVVEDLLRAGADLSLLDRLGNSVLHLAAKEGHDKVLSILLKHKKAALLLDHPNGDGLNAIHLAMMSNSLPCLLLLVAAGADVNAQEQKSGRTALHLAVEHDNISLAGCLLLEGDAHVDSTTYDGTTPLHIAAGRGSTRLAALLKAAGADPLVENFEPLYDLDDSWENAGEDEGVVPGTTPLDMATSWQVFDILNGKPYEPEFTSDDLLAQGDMKQLAEDVKLQLYKLLEIPDPDKNWATLAQKLGLGILNNAFRLSPAPSKTLMDNYEVSGGTVRELVEALRQMGYTEAIEVIQAASSPVKTTSQAHSLPLSPASTRQQIDELRDSDSVCDSGVETSFRKLSFTESLTSGASLLTLNKMPHDYGQEGPLEGKIXXMAEDDPYLGRPEQMFHLDPSLTHTIFNPEVFQPQMALPTADGPYLQILEQPKQRGFRFRYVCEGPSHGGLPGASSEKNKKSYPQVKICNYVGPAKVIVQLVTNGKNIHLHAHSLVGKHCEDGICTVTAGPKDMVVGFANLGILHVTKKKVFETLEARMTEACIRGYNPGLLVHPDLAYLQAEGGGDRQLGDREKELIRQAALQQTKEMDLSVVRLMFTAFLPDSTGSFTRRLEPVVSDAIYDSKAPNASNLKIVRMDRTAGCVTGGEEIYLLCDKVQKDDIQIRFYEEEENGGVWEGFGDFSPTDVHRQFAIVFKTPKYKDINITKPASVFVQLRRKSDLETSEPKPFLYYPEIKDKEEVQRKRQKLMPNFSDSFGGGSGAGAGGGGMFGSGGGGGGTGSTGPGYSFPHYGFPTYGGITFHPGTTKSNAGMKHGTMDTESKKDPEGCDKSDDKNTVNLFGKVIETTEQDQEPSEATVGNGEVTLTYATGTKEESAGVQDNLFLEKAMQLAKRHANALFDYAVTGDVKMLLAVQRHLTAVQDENGDSVLHLAIIHLHSQLVRDLLEVTSGLISDDIINMRNDLYQTPLHLAVITKQEDVVEDLLRAGADLSLLDRLGNSVLHLAAKEGHDKVLSILLKHKKAALLLDHPNGDGLNAIHLAMMSNSLPCLLLLVAAGADVNAQEQKSGRTALHLAVEHDNISLAGCLLLEGDAHVDSTTYDGTTPLHIAAGRGSTRLAALLKAAGADPLVENFEPLYDLDDSWENAGEDEGVVPGTTPLDMATSWQVFDILNGKPYEPEFTSDDLLAQGDMKQLAEDVKLQLYKLLEIPDPDKNWATLAQKLGLGILNNAFRLSPAPSKTLMDNYEVSGGTVRELVEALRQMGYTEAIEVIQAASSPVKTTSQAHSLPLSPASTRQQIDELRDSDSVCDSGVETSFRKLSFTESLTSGASLLTLNKMPHDYGQEGPLEGKIXXMPRCPAGAMDEGPVDLRTRPKAAGLPGAALPLRKRPLRAPSPEPAAPRGAAGLVVPLDPLRGGCDLPAVPGPPHGLARPEALYYPGALLPLYPTRAMGSPFPLVNLPTPLYPMMCPMEHPLSADIAMATRADEDGDTPLHIAVVQGNLPAVHRLVNLFQQGGRELDIYNNLRQTPLHLAVITTLPSVVRLLVTAGASPMALDRHGQTAAHLACEHRSPTCLRALLDSAAPGTLDLEARNYDGLTALHVAVNTECQETVQLLLERGADIDAVDIKSGRSPLIHAVENNSLSMVQLLLQHGANVNAQMYSGSSALHSASGRGLLPLVRTLVRSGADSSLKNCHNDTPLMVARSRRVIDILRGKATRPASTSQPDPSPDRSANTSPESSSRLSSNGLLSASPSSSPSQSPPRDPPGFPMAPPNFFLPSPSPPAFLPFAGVLRGPGRPVPPSPAPGGSXXMAEDDPYLGRPEQMFHLDPSLTHTIFNPEVFQPQMALPTDGPYLQILEQPKQRGFRFRYVCEGPSHGGLPGASSEKNKKSYPQVKICNYVGPAKVIVQLVTNGKNIHLHAHSLVGKHCEDGICTVTAGPKDMVVGFANLGILHVTKKKVFETLEARMTEACIRGYNPGLLVHPDLAYLQAEGGGDRQLGDREKELIRQAALQQTKEMDLSVVRLMFTAFLPDSTGSFTRRLEPVVSDAIYDSKAPNASNLKIVRMDRTAGCVTGGEEIYLLCDKVQKDDIQIRFYEEEENGGVWEGFGDFSPTDVHRQFAIVFKTPKYKDINITKPASVFVQLRRKSDLETSEPKPFLYYPEIKDKEEVQRKRQKLMPNFSDSFGGGSGAGAGGGGMFGSGGGGGGTGSTGPGYSFPHYGFPTYGGITFHPGTTKSNAGMKHGTMDTESKKDPEGCDKSDDKNTVNLFGKVIETTEQDQEPSEATVGNGEVTLTYATGTKEESAGVQDNLFLEKAMQLAKRHANALFDYAVTGDVKMLLAVQRHLTAVQDENGDSVLHLAIIHLHSQLVRDLLEVTSGLISDDIINMRNDLYQTPLHLAVITKQEDVVEDLLRAGADLSLLDRLGNSVLHLAAKEGHDKVLSILLKHKKAALLLDHPNGDGLNAIHLAMMSNSLPCLLLLVAAGADVNAQEQKSGRTALHLAVEHDNISLAGCLLLEGDAHVDSTTYDGTTPLHIAAGRGSTRLAALLKAAGADPLVENFEPLYDLDDSWENAGEDEGVVPGTTPLDMATSWQVFDILNGKPYEPEFTSDDLLAQGDMKQLAEDVKLQLYKLLEIPDPDKNWATLAQKLGLGILNNAFRLSPAPSKTLMDNYEVSGGTVRELVEALRQMGYTEAIEVIQAASSPVKTTSQAHSLPLSPASTRQQIDELRDSDSVCDSGVETSFRKLSFTESLTSGASLLTLNKMPHDYGQEGPLEGKIX"
##      Mutant_Peptide         Total_RNA Tumor_RNA_Ratio Tumor_RNA
## [1,] "MAEDDPYLGRPEQMFHLILL" "NA"      "NA"            "NA"     
## [2,] "MAEDDPYLGRPEQMFHLILL" "NA"      "NA"            "NA"     
## [3,] "MAEDDPYLGRPEQMFHLILL" "NA"      "NA"            "NA"     
## [4,] "MAEDDPYLGRPEQMFHLILL" "NA"      "NA"            "NA"     
## [5,] "MAEDDPYLGRPEQMFHLILL" "NA"      "NA"            "NA"     
## [6,] "MAEDDPYLGRPEQMFHLILL" "NA"      "NA"            "NA"     
##      Tumor_RNA_based_on_DNA MutRatio MutRatio_Min MutRatio_Max
## [1,] "NA"                   "NA"     "NA"         "NA"        
## [2,] "NA"                   "NA"     "NA"         "NA"        
## [3,] "NA"                   "NA"     "NA"         "NA"        
## [4,] "NA"                   "NA"     "NA"         "NA"        
## [5,] "NA"                   "NA"     "NA"         "NA"        
## [6,] "NA"                   "NA"     "NA"         "NA"
```

```
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
## 
## [[2]]
##                                                        
## Num_Peptide_Per_NM                                3.000
## Num_Cond_Peptide_Per_NM                           3.000
## Num_Rest_Peptide_Per_NM                           1.000
## Num_Rest_Peptide_Per_NM / Num_Cond_Peptide_Per_NM 0.333
## -logP                                             0.460
##
## [[3]]
##                                                     MAEDDPYLGRPEQMFHLILL-0_atggcagaag-0_1
## Num_Peptide_Per_Pep                                                                 3.000
## Num_Cond_Peptide_Per_Pep                                                            3.000
## Num_Rest_Peptide_Per_Pep                                                            1.000
## Num_Rest_Peptide_Per_Pep / Num_Cond_Peptide_Per_Pep                                 0.333
## -logP                                                                               0.460
```