---
title: "SampleCodeForNeoantimon"
author: "T. Hasegawa"
date: "2018/05/02"
output: html_document
---

<!--
%\VignetteEngine{knitr::knitr}
%\VignetteIndexEntry{Sample Code to Use Neoantimon}
-->

## Data Preparation and Sample Codes for Analysis

```r
#install.packages('devtools');
library(devtools);
install_github('hase62/Neoantimon');
```

```
## Downloading GitHub repo hase62/Neoantimon@master
## from URL https://api.github.com/repos/hase62/Neoantimon/zipball/master
```

```
## Installing Neoantimon
```

```
## '/Library/Frameworks/R.framework/Resources/bin/R' --no-site-file  \
##   --no-environ --no-save --no-restore --quiet CMD INSTALL  \
##   '/private/var/folders/s3/nbtvmrx149v02tzgjnm2_6nc0000gn/T/RtmpTyR9zp/devtools23d24a5ce764/hase62-Neoantimon-a8454c9'  \
##   --library='/Library/Frameworks/R.framework/Versions/3.4/Resources/library'  \
##   --install-tests
```

```
## 
```

```r
library(Neoantimon);
```


```r
data("sample_vcf")
head(sample_vcf, row.names = FALSE)
```

```
##   Chr     Start       End Ref Alt Func.refGene Gene.refGene
## 1   1  47399872  47399872   A   G       exonic      CYP4A11
## 2   1 116941338 116941338   T   C       exonic       ATP1A1
## 3   4  24556416  24556416   T   C       exonic        DHX15
## 4   4  70156404  70156404   -   T       exonic      UGT2B28
## 5   6  75899298  75899298   T   -       exonic      COL12A1
## 6   9  89561162  89561162   C   T       exonic         GAS1
##   GeneDetail.refGene ExonicFunc.refGene
## 1      nonsynonymous                SNV
## 2         synonymous                SNV
## 3      nonsynonymous                SNV
## 4         frameshift          insertion
## 5         frameshift           deletion
## 6      nonsynonymous                SNV
##                                                                                                                   AAChange.refGene
## 1                                                                                         CYP4A11:NM_000778:exon8:c.T1064C:p.L355P
## 2 ATP1A1:NM_000701:exon16:c.T2220C:p.D740D,ATP1A1:NM_001160233:exon16:c.T2220C:p.D740D,ATP1A1:NM_001160234:exon16:c.T2127C:p.D709D
## 3                                                                                           DHX15:NM_001358:exon5:c.A1012G:p.T338A
## 4                                                                                      UGT2B28:NM_053039:exon5:c.1186dupT:p.L395fs
## 5                                                                                       COL12A1:NM_004370:exon6:c.628delA:p.I210fs
## 6                                                                                             GAS1:NM_002048:exon1:c.G533A:p.R178H
##   cytoBand depth_tumor variantNum_tumor depth_normal
## 1     1p33          64               28           49
## 2   1p13.1         100               39          111
## 3   4p15.2         143               47          151
## 4   4q13.2          43               15           41
## 5     6q13         122               38           73
## 6  9q21.33          20                5           26
```


```r
data("sample_sv_bnd")
head(sample_sv_bnd, row.names = FALSE)
```

```
##   Chr     Start       End Ref            Alt Func.refGene Gene.refGene
## 1   1 115005805 115005805   C C]20:34827929]       exonic       TRIM33
## 2   1 204908711 204908711   A [2:172743385[A     intronic        NFASC
## 3   2  25870534  25870534   T  ]2:25965720]G     intronic         DTNB
## 4   2  25965720  25965720   T  T[2:25870536[       exonic        ASXL2
## 5   2 214794791 214794791   C C[2:214798169[       exonic       SPAG16
## 6   2 214798169 214798169   T ]2:214794791]T     intronic       SPAG16
##         mateID
## 1 SVMERGE137_1
## 2  SVMERGE15_1
## 3 SVMERGE116_1
## 4 SVMERGE116_2
## 5   SVMERGE3_1
## 6   SVMERGE3_2
```


```r
library(Neoantimon)
data("sample_hla_table_c1")
head(sample_hla_table_c1, row.names = FALSE)
```

```
##      Name      A1      A2      B1      B2      C1      C2
## 1  sample A*02:01 A*32:01 B*15:17 B*51:01 C*07:01 C*15:02
## 2 sample2 A*02:01 A*32:01 B*15:17 B*51:01 C*07:01 C*15:02
```


```r
data("sample_hla_table_c2")
head(sample_hla_table_c2, row.names = FALSE)
```

```
##      Name      DPA11      DPA12      DPB11      DPB12      DQA11
## 1  sample DPA1*01:03 DPA1*02:01 DPB1*02:01 DPB1*09:01 DQA1*01:02
## 2 sample2 DPA1*01:03 DPA1*02:01 DPB1*02:01 DPB1*09:01 DQA1*01:02
##        DQA12      DQB11      DQB12      DRB11      DRB12
## 1 DQA1*05:05 DQB1*03:01 DQB1*06:04 DRB1*11:04 DRB1*13:02
## 2 DQA1*05:05 DQB1*03:01 DQB1*06:04 DRB1*11:04 DRB1*13:02
```


```r
data("sample_rna_exp")
head(sample_rna_exp, row.names = FALSE)
```

```
##   gene_short_name                    locus expression
## 1         5S_rRNA     12:34358633-34358752    0.89094
## 2         5S_rRNA     14:68123821-68123928    8.34293
## 3         5S_rRNA     16:34977638-34990886    6.47095
## 4         5S_rRNA    2:162266064-162266181    0.00000
## 5         5S_rRNA GL000192.1:415331-415454    1.04608
## 6         5S_rRNA   GL000228.1:20112-20230    0.22506
```


```r
data("sample_copynum")
head(sample_copynum, row.names = FALSE)
```

```
##               Chromosome Position      Log.R segmented.LogR BAF
## SNP_A-8575125          1   564621  0.6071447    -0.09862298   1
## SNP_A-2205441          1   799463  0.1519967    -0.09862298   1
## SNP_A-8325638          1  1017216  0.8146911    -0.09862298   0
## SNP_A-8450345          1  1158277 -1.9594627    -0.09862298   0
## SNP_A-8507155          1  1242215  0.2927962    -0.09862298   0
## SNP_A-8596426          1  1462766 -0.2234090    -0.09862298   1
##               segmented.BAF Copy.number Minor.allele Raw.copy.number
## SNP_A-8575125            NA           2            1       4.3540752
## SNP_A-2205441            NA           2            1       2.1467339
## SNP_A-8325638            NA           3            1       5.8658499
## SNP_A-8450345            NA           1            0      -0.5035897
## SNP_A-8507155            NA           2            1       2.6999875
## SNP_A-8596426            NA           2            1       1.0726687
```


```r
data("sample_result_SNV_CLASS1_ALL")
head(sample_result_SNV_CLASS1_ALL, row.names = FALSE)
```

```
##           HLA Pos      Gene Evaluated_Mutant_Peptide Mut_IC50 Mut_Rank
## 1 HLA-A*02:01   2 0_CYP4A11            HQERCREEIHSLP  37362.9       55
## 2 HLA-A*02:01   3 0_CYP4A11            QERCREEIHSLPG  45831.8       90
## 3 HLA-A*02:01   4 0_CYP4A11            ERCREEIHSLPGD  45989.7       90
## 4 HLA-A*02:01   5 0_CYP4A11            RCREEIHSLPGDG  46482.0       95
## 5 HLA-A*02:01   6 0_CYP4A11            CREEIHSLPGDGA  42936.1       80
## 6 HLA-A*02:01   7 0_CYP4A11            REEIHSLPGDGAS  44359.4       85
##   Evaluated_Wt_Peptide Wt_IC50 Wt_Rank Chr     NM_ID   Change Ref Alt Prob
## 1        HQERCREEIHSLL 27284.3      32   1 NM_000778 c.T1064C   A   G    0
## 2        QERCREEIHSLLG 44161.6      85   1 NM_000778 c.T1064C   A   G    0
## 3        ERCREEIHSLLGD 44575.4      85   1 NM_000778 c.T1064C   A   G    0
## 4        RCREEIHSLLGDG 44377.6      85   1 NM_000778 c.T1064C   A   G    0
## 5        CREEIHSLLGDGA 41760.9      70   1 NM_000778 c.T1064C   A   G    0
## 6        REEIHSLLGDGAS 43379.7      80   1 NM_000778 c.T1064C   A   G    0
##   Mutation_Prob. Exon_Start Exon_End Mutation_Position Total_Depth
## 1              0   47394845 47407156        1_47399872         113
## 2              0   47394845 47407156        1_47399872         113
## 3              0   47394845 47407156        1_47399872         113
## 4              0   47394845 47407156        1_47399872         113
## 5              0   47394845 47407156        1_47399872         113
## 6              0   47394845 47407156        1_47399872         113
##   Tumor_Depth                      Wt_Peptide
## 1          64 HPKHQERCREEIHSLLGDGASITWNHLDQMP
## 2          64 HPKHQERCREEIHSLLGDGASITWNHLDQMP
## 3          64 HPKHQERCREEIHSLLGDGASITWNHLDQMP
## 4          64 HPKHQERCREEIHSLLGDGASITWNHLDQMP
## 5          64 HPKHQERCREEIHSLLGDGASITWNHLDQMP
## 6          64 HPKHQERCREEIHSLLGDGASITWNHLDQMP
##                    Mutant_Peptide Total_RNA Tumor_RNA_Ratio Tumor_RNA
## 1 HPKHQERCREEIHSLPGDGASITWNHLDQMP        NA              NA        NA
## 2 HPKHQERCREEIHSLPGDGASITWNHLDQMP        NA              NA        NA
## 3 HPKHQERCREEIHSLPGDGASITWNHLDQMP        NA              NA        NA
## 4 HPKHQERCREEIHSLPGDGASITWNHLDQMP        NA              NA        NA
## 5 HPKHQERCREEIHSLPGDGASITWNHLDQMP        NA              NA        NA
## 6 HPKHQERCREEIHSLPGDGASITWNHLDQMP        NA              NA        NA
##   Tumor_RNA_based_on_DNA MutRatio MutRatio_Min MutRatio_Max
## 1                     NA       NA           NA           NA
## 2                     NA       NA           NA           NA
## 3                     NA       NA           NA           NA
## 4                     NA       NA           NA           NA
## 5                     NA       NA           NA           NA
## 6                     NA       NA           NA           NA
```


```r
data("sample_result_SNV_CLASS2_ALL")
head(sample_result_SNV_CLASS2_ALL, row.names = FALSE)
```

```
##                     HLA Pos      Gene Evaluated_Mutant_Peptide Mut_IC50
## 1 HLA-DPA10103-DPB10201   1 0_CYP4A11          PKHQERCREEIHSLP  9714.36
## 2 HLA-DPA10103-DPB10201   2 0_CYP4A11          KHQERCREEIHSLPG  9186.08
## 3 HLA-DPA10103-DPB10201   3 0_CYP4A11          HQERCREEIHSLPGD  9621.48
## 4 HLA-DPA10103-DPB10201   4 0_CYP4A11          QERCREEIHSLPGDG 10479.94
## 5 HLA-DPA10103-DPB10201   5 0_CYP4A11          ERCREEIHSLPGDGA 10729.49
## 6 HLA-DPA10103-DPB10201   6 0_CYP4A11          RCREEIHSLPGDGAS 11385.65
##   Mut_Rank Evaluated_Wt_Peptide Wt_IC50 Wt_Rank Chr     NM_ID   Change Ref
## 1       90      PKHQERCREEIHSLL 2610.11      60   1 NM_000778 c.T1064C   A
## 2       90      KHQERCREEIHSLLG 2097.09      55   1 NM_000778 c.T1064C   A
## 3       90      HQERCREEIHSLLGD 1992.73      55   1 NM_000778 c.T1064C   A
## 4       90      QERCREEIHSLLGDG 2154.36      55   1 NM_000778 c.T1064C   A
## 5       90      ERCREEIHSLLGDGA 2193.35      60   1 NM_000778 c.T1064C   A
## 6       95      RCREEIHSLLGDGAS 2362.58      60   1 NM_000778 c.T1064C   A
##   Alt Prob Mutation_Prob. Exon_Start Exon_End Mutation_Position
## 1   G    0              0   47394845 47407156        1_47399872
## 2   G    0              0   47394845 47407156        1_47399872
## 3   G    0              0   47394845 47407156        1_47399872
## 4   G    0              0   47394845 47407156        1_47399872
## 5   G    0              0   47394845 47407156        1_47399872
## 6   G    0              0   47394845 47407156        1_47399872
##   Total_Depth Tumor_Depth                      Wt_Peptide
## 1         113          64 HPKHQERCREEIHSLLGDGASITWNHLDQMP
## 2         113          64 HPKHQERCREEIHSLLGDGASITWNHLDQMP
## 3         113          64 HPKHQERCREEIHSLLGDGASITWNHLDQMP
## 4         113          64 HPKHQERCREEIHSLLGDGASITWNHLDQMP
## 5         113          64 HPKHQERCREEIHSLLGDGASITWNHLDQMP
## 6         113          64 HPKHQERCREEIHSLLGDGASITWNHLDQMP
##                    Mutant_Peptide Total_RNA Tumor_RNA_Ratio Tumor_RNA
## 1 HPKHQERCREEIHSLPGDGASITWNHLDQMP        NA              NA        NA
## 2 HPKHQERCREEIHSLPGDGASITWNHLDQMP        NA              NA        NA
## 3 HPKHQERCREEIHSLPGDGASITWNHLDQMP        NA              NA        NA
## 4 HPKHQERCREEIHSLPGDGASITWNHLDQMP        NA              NA        NA
## 5 HPKHQERCREEIHSLPGDGASITWNHLDQMP        NA              NA        NA
## 6 HPKHQERCREEIHSLPGDGASITWNHLDQMP        NA              NA        NA
##   Tumor_RNA_based_on_DNA MutRatio MutRatio_Min MutRatio_Max
## 1                     NA       NA           NA           NA
## 2                     NA       NA           NA           NA
## 3                     NA       NA           NA           NA
## 4                     NA       NA           NA           NA
## 5                     NA       NA           NA           NA
## 6                     NA       NA           NA           NA
```


```r
data("sample_result_INDEL_CLASS1_ALL")
head(sample_result_INDEL_CLASS1_ALL, row.names = FALSE)
```

```
##           HLA Pos      Gene Evaluated_Mutant_Peptide_Core
## 1 HLA-A*02:01   1 0_UGT2B28                    GIPMVGIPLV
## 2 HLA-A*02:01   2 0_UGT2B28                 YHGIPMVGIPLVL
## 3 HLA-A*02:01   3 0_UGT2B28                  HGIPMVGIPLVL
## 4 HLA-A*02:01   4 0_UGT2B28                    GIPMVGIPLV
## 5 HLA-A*02:01   5 0_UGT2B28                    IPMVGIPLVL
## 6 HLA-A*02:01   2 0_UGT2B28                    GIPMVGIPLV
##   Evaluated_Mutant_Peptide Mut_IC50 Mut_Rank Chr     NM_ID         Change
## 1            IYHGIPMVGIPLV   4359.5      8.0   4 NM_053039 Out_c.1186dupT
## 2            YHGIPMVGIPLVL    482.1      3.0   4 NM_053039 Out_c.1186dupT
## 3            HGIPMVGIPLVLG  10645.4     14.0   4 NM_053039 Out_c.1186dupT
## 4            GIPMVGIPLVLGS  10291.2     14.0   4 NM_053039 Out_c.1186dupT
## 5            IPMVGIPLVLGST  13857.7     17.0   4 NM_053039 Out_c.1186dupT
## 6             YHGIPMVGIPLV   3641.2      7.5   4 NM_053039 Out_c.1186dupT
##   Ref Alt Prob Mutation_Prob. Exon_Start Exon_End Mutation_Position
## 1   -   T    0              0   70146216 70160768        4_70156404
## 2   -   T    0              0   70146216 70160768        4_70156404
## 3   -   T    0              0   70146216 70160768        4_70156404
## 4   -   T    0              0   70146216 70160768        4_70156404
## 5   -   T    0              0   70146216 70160768        4_70156404
## 6   -   T    0              0   70146216 70160768        4_70156404
##   Total_Depth Tumor_Depth                     Wt_Peptide
## 1          84          43 EAIYHGIPMVGIPLFWDQPDNIAHMKAKGA
## 2          84          43 EAIYHGIPMVGIPLFWDQPDNIAHMKAKGA
## 3          84          43 EAIYHGIPMVGIPLFWDQPDNIAHMKAKGA
## 4          84          43 EAIYHGIPMVGIPLFWDQPDNIAHMKAKGA
## 5          84          43 EAIYHGIPMVGIPLFWDQPDNIAHMKAKGA
## 6          84          43 EAIYHGIPMVGIPLFWDQPDNIAHMKAKGA
##         Mutant_Peptide Total_RNA Tumor_RNA_Ratio Tumor_RNA
## 1 EAIYHGIPMVGIPLVLGSTX        NA              NA        NA
## 2 EAIYHGIPMVGIPLVLGSTX        NA              NA        NA
## 3 EAIYHGIPMVGIPLVLGSTX        NA              NA        NA
## 4 EAIYHGIPMVGIPLVLGSTX        NA              NA        NA
## 5 EAIYHGIPMVGIPLVLGSTX        NA              NA        NA
## 6 EAIYHGIPMVGIPLVLGSTX        NA              NA        NA
##   Tumor_RNA_based_on_DNA MutRatio MutRatio_Min MutRatio_Max
## 1                     NA       NA           NA           NA
## 2                     NA       NA           NA           NA
## 3                     NA       NA           NA           NA
## 4                     NA       NA           NA           NA
## 5                     NA       NA           NA           NA
## 6                     NA       NA           NA           NA
```


```r
data("sample_result_INDEL_CLASS2_ALL")
head(sample_result_INDEL_CLASS2_ALL)
```

```
##                     HLA Pos      Gene Evaluated_Mutant_Peptide_Core
## 1 HLA-DPA10103-DPB10201   0 0_UGT2B28                     YHGIPMVGI
## 2 HLA-DPA10103-DPB10201   1 0_UGT2B28                     YHGIPMVGI
## 3 HLA-DPA10103-DPB10201   2 0_UGT2B28                     MVGIPLVLG
## 4 HLA-DPA10103-DPB10201   3 0_UGT2B28                     MVGIPLVLG
## 5 HLA-DPA10103-DPB10201   4 0_UGT2B28                     MVGIPLVLG
## 6 HLA-DPA10103-DPB10201   0 1_COL12A1                     YQRDELLAA
##   Evaluated_Mutant_Peptide Mut_IC50 Mut_Rank Chr     NM_ID         Change
## 1          EAIYHGIPMVGIPLV  1028.70     40.0   4 NM_053039 Out_c.1186dupT
## 2          AIYHGIPMVGIPLVL  1018.07     40.0   4 NM_053039 Out_c.1186dupT
## 3          IYHGIPMVGIPLVLG  1368.82     46.0   4 NM_053039 Out_c.1186dupT
## 4          YHGIPMVGIPLVLGS  1613.59     49.0   4 NM_053039 Out_c.1186dupT
## 5          HGIPMVGIPLVLGST  2001.80     55.0   4 NM_053039 Out_c.1186dupT
## 6          QYYQRDELLAAIKKF    93.99      8.5   6 NM_004370  Out_c.628delA
##   Ref Alt Prob Mutation_Prob. Exon_Start Exon_End Mutation_Position
## 1   -   T    0              0   70146216 70160768        4_70156404
## 2   -   T    0              0   70146216 70160768        4_70156404
## 3   -   T    0              0   70146216 70160768        4_70156404
## 4   -   T    0              0   70146216 70160768        4_70156404
## 5   -   T    0              0   70146216 70160768        4_70156404
## 6   T   -    0              0   75794041 75915623        6_75899298
##   Total_Depth Tumor_Depth                        Wt_Peptide
## 1          84          43    EAIYHGIPMVGIPLFWDQPDNIAHMKAKGA
## 2          84          43    EAIYHGIPMVGIPLFWDQPDNIAHMKAKGA
## 3          84          43    EAIYHGIPMVGIPLFWDQPDNIAHMKAKGA
## 4          84          43    EAIYHGIPMVGIPLFWDQPDNIAHMKAKGA
## 5          84          43    EAIYHGIPMVGIPLFWDQPDNIAHMKAKGA
## 6         195         122 QYYQRDELLAAIKKIPYKGGNTMTGDAIDYLVK
##            Mutant_Peptide Total_RNA Tumor_RNA_Ratio Tumor_RNA
## 1    EAIYHGIPMVGIPLVLGSTX        NA              NA        NA
## 2    EAIYHGIPMVGIPLVLGSTX        NA              NA        NA
## 3    EAIYHGIPMVGIPLVLGSTX        NA              NA        NA
## 4    EAIYHGIPMVGIPLVLGSTX        NA              NA        NA
## 5    EAIYHGIPMVGIPLVLGSTX        NA              NA        NA
## 6 QYYQRDELLAAIKKFHIKVATQX        NA              NA        NA
##   Tumor_RNA_based_on_DNA MutRatio MutRatio_Min MutRatio_Max
## 1                     NA       NA           NA           NA
## 2                     NA       NA           NA           NA
## 3                     NA       NA           NA           NA
## 4                     NA       NA           NA           NA
## 5                     NA       NA           NA           NA
## 6                     NA       NA           NA           NA
```










