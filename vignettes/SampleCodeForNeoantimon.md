---
title: "SampleCodeForNeoantimon"
author: "T. Hasegawa"
date: "2024/6/12"
output:
pdf_document: default
html_document: default
---

```{=html}
<!--
%\VignetteEngine{knitr::knitr}
%\VignetteIndexEntry{Sample Code to Use Neoantimon}
-->
```
## Data Preparation and Sample Codes for Analysis


``` r
#install.packages('devtools');
library(devtools);
```

```
## Loading required package: usethis
```

``` r
install_github('hase62/Neoantimon');
```

```
## Skipping install of 'Neoantimon' from a github remote, the SHA1 (2ed545f8) has not changed since last install.
##   Use `force = TRUE` to force installation
```

``` r
library(Neoantimon);
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("biomaRt")
```

```
## Bioconductor version 3.19 (BiocManager 1.30.23), R 4.4.0 (2024-04-24)
```

```
## Warning: package(s) not installed when version(s) same as or greater than current; use
##   `force = TRUE` to re-install: 'biomaRt'
```

``` r
library(biomaRt)
```

To calculate the binding affinity of neoantigen candaites, which are generated from from SNVs, to HLA ClassI. When using MHCflurry, [[1]] and [[2]] include the results of NetMHCpan and MHCflurry, respectively.


``` r
Result_HLA1_SNV <- MainSNVClass1(input_annovar_format_file = "./../data/sample_vcf.annovar.txt",
                                   file_name_in_hla_table = "sample",
                                   hla_file = "./../data/sample_hla_table_c1.txt",
                                   refflat_file  = "./../lib/refFlat.grch37.txt",
                                   refmrna_file = "./../lib/refMrna.grch37.fa",
                                   rnaexp_file = "./../data/sample_rna_exp.txt",
                                   netMHCpan_dir = "./../lib/netMHCpan-4.1/netMHCpan",
                                   depth_tumor_column = 12,
                                   depth_normal_column = 14,
                                   SNPs = "./../data/sample_vcf.snps.txt",
                                   multiple_variants = TRUE)
```

```
## [1] "HLAtype: A*02:01" "HLAtype: A*32:01" "HLAtype: B*15:17" "HLAtype: B*51:01"
## [5] "HLAtype: C*07:01" "HLAtype: C*15:02"
## [1] "Set chr_column as 1 Chr"
## [1] "Set mutation_start_column as 2 Start"
## [1] "Set mutation_end_column as 3 End"
## [1] "Set mutation_ref_column as 4 Ref"
## [1] "Set mutation_alt_column as 5 Alt"
## [1] "Set depth_normal_column as 14 depth_normal"
## [1] "Set depth_tumor_column as 12 depth_tumor"
## [1] "Set nm_id_column as 10 AAChange.refGene"
## [1] "Please Confirm that Reading Start Line is 1"
## [1] "Chr\tStart\tEnd\tRef\tAlt\tFunc.refGene\tGene.refGene\tGeneDetail.refGene\tExonicFunc.refGene\tAAChange.refGene\tcytoBand\tdepth_tumor\tvariantNum_tumor\tdepth_normal"
## [1] "Please Confirm that Reading Start Line is 1"
## [1] "CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT"
## [1] "Start Analysis: Mutation 1"
## [1] "NM_001358 - Ref Matched to and NM_ID Attached Information"
## [1] "Peptide Successfully Generated!!"
## [1] "Start Analysis: Mutation 2"
## [1] "NM_002048 - Ref Matched to and NM_ID Attached Information"
## [1] "Peptide Successfully Generated!!"
## [1] "Start Analysis: Mutation 3"
## [1] "NM_006205 - Ref Matched to and NM_ID Attached Information"
## [1] "Peptide Successfully Generated!!"
## [1] "result.ID.SNV1/sample_vcf.annovar.txt.ID_SNV.peptide.txt.list.vcf"
## [1] "Executing netMHCpan to result.ID.SNV1"
## [1] "./../lib/netMHCpan-4.1/netMHCpan -BA  -l 8,9,10,11 -f result.ID.SNV1/sample_vcf.annovar.txt.ID_SNV.peptide.fasta -a HLA-A02:01 > result.ID.SNV1/sample_vcf.annovar.txt.ID_SNV.HLACLASS1.1.peptide.txt"
## [1] "./../lib/netMHCpan-4.1/netMHCpan -BA  -l 8,9,10,11 -f result.ID.SNV1/sample_vcf.annovar.txt.ID_SNV.peptide.fasta -a HLA-A32:01 > result.ID.SNV1/sample_vcf.annovar.txt.ID_SNV.HLACLASS1.2.peptide.txt"
## [1] "./../lib/netMHCpan-4.1/netMHCpan -BA  -l 8,9,10,11 -f result.ID.SNV1/sample_vcf.annovar.txt.ID_SNV.peptide.fasta -a HLA-B15:17 > result.ID.SNV1/sample_vcf.annovar.txt.ID_SNV.HLACLASS1.3.peptide.txt"
## [1] "./../lib/netMHCpan-4.1/netMHCpan -BA  -l 8,9,10,11 -f result.ID.SNV1/sample_vcf.annovar.txt.ID_SNV.peptide.fasta -a HLA-B51:01 > result.ID.SNV1/sample_vcf.annovar.txt.ID_SNV.HLACLASS1.4.peptide.txt"
## [1] "./../lib/netMHCpan-4.1/netMHCpan -BA  -l 8,9,10,11 -f result.ID.SNV1/sample_vcf.annovar.txt.ID_SNV.peptide.fasta -a HLA-C07:01 > result.ID.SNV1/sample_vcf.annovar.txt.ID_SNV.HLACLASS1.5.peptide.txt"
## [1] "./../lib/netMHCpan-4.1/netMHCpan -BA  -l 8,9,10,11 -f result.ID.SNV1/sample_vcf.annovar.txt.ID_SNV.peptide.fasta -a HLA-C15:02 > result.ID.SNV1/sample_vcf.annovar.txt.ID_SNV.HLACLASS1.6.peptide.txt"
## [1] "./../lib/netMHCpan-4.1/netMHCpan -BA  -l 8,9,10,11 -f result.ID.SNV1/sample_vcf.annovar.txt.ID_SNV.wtpeptide.fasta -a HLA-A02:01 > result.ID.SNV1/sample_vcf.annovar.txt.ID_SNV.HLACLASS1.1.wtpeptide.txt"
## [1] "./../lib/netMHCpan-4.1/netMHCpan -BA  -l 8,9,10,11 -f result.ID.SNV1/sample_vcf.annovar.txt.ID_SNV.wtpeptide.fasta -a HLA-A32:01 > result.ID.SNV1/sample_vcf.annovar.txt.ID_SNV.HLACLASS1.2.wtpeptide.txt"
## [1] "./../lib/netMHCpan-4.1/netMHCpan -BA  -l 8,9,10,11 -f result.ID.SNV1/sample_vcf.annovar.txt.ID_SNV.wtpeptide.fasta -a HLA-B15:17 > result.ID.SNV1/sample_vcf.annovar.txt.ID_SNV.HLACLASS1.3.wtpeptide.txt"
## [1] "./../lib/netMHCpan-4.1/netMHCpan -BA  -l 8,9,10,11 -f result.ID.SNV1/sample_vcf.annovar.txt.ID_SNV.wtpeptide.fasta -a HLA-B51:01 > result.ID.SNV1/sample_vcf.annovar.txt.ID_SNV.HLACLASS1.4.wtpeptide.txt"
## [1] "./../lib/netMHCpan-4.1/netMHCpan -BA  -l 8,9,10,11 -f result.ID.SNV1/sample_vcf.annovar.txt.ID_SNV.wtpeptide.fasta -a HLA-C07:01 > result.ID.SNV1/sample_vcf.annovar.txt.ID_SNV.HLACLASS1.5.wtpeptide.txt"
## [1] "./../lib/netMHCpan-4.1/netMHCpan -BA  -l 8,9,10,11 -f result.ID.SNV1/sample_vcf.annovar.txt.ID_SNV.wtpeptide.fasta -a HLA-C15:02 > result.ID.SNV1/sample_vcf.annovar.txt.ID_SNV.HLACLASS1.6.wtpeptide.txt"
## [1] "Merging Results..."
## [1] "/Users/takanorihasegawa/Documents/GitHub/Neoantimon/vignettes/result.ID.SNV1/sample_vcf.annovar.txt.ID_SNV.HLACLASS1.1.peptide.txt"
## [1] "33.3333333333333 perc. fin"
## [1] "66.6666666666667 perc. fin"
## [1] "100 perc. fin"
## [1] "/Users/takanorihasegawa/Documents/GitHub/Neoantimon/vignettes/result.ID.SNV1/sample_vcf.annovar.txt.ID_SNV.HLACLASS1.2.peptide.txt"
## [1] "33.3333333333333 perc. fin"
## [1] "66.6666666666667 perc. fin"
## [1] "100 perc. fin"
## [1] "/Users/takanorihasegawa/Documents/GitHub/Neoantimon/vignettes/result.ID.SNV1/sample_vcf.annovar.txt.ID_SNV.HLACLASS1.3.peptide.txt"
## [1] "33.3333333333333 perc. fin"
## [1] "66.6666666666667 perc. fin"
## [1] "100 perc. fin"
## [1] "/Users/takanorihasegawa/Documents/GitHub/Neoantimon/vignettes/result.ID.SNV1/sample_vcf.annovar.txt.ID_SNV.HLACLASS1.4.peptide.txt"
## [1] "33.3333333333333 perc. fin"
## [1] "66.6666666666667 perc. fin"
## [1] "100 perc. fin"
## [1] "/Users/takanorihasegawa/Documents/GitHub/Neoantimon/vignettes/result.ID.SNV1/sample_vcf.annovar.txt.ID_SNV.HLACLASS1.5.peptide.txt"
## [1] "33.3333333333333 perc. fin"
## [1] "66.6666666666667 perc. fin"
## [1] "100 perc. fin"
## [1] "/Users/takanorihasegawa/Documents/GitHub/Neoantimon/vignettes/result.ID.SNV1/sample_vcf.annovar.txt.ID_SNV.HLACLASS1.6.peptide.txt"
## [1] "33.3333333333333 perc. fin"
## [1] "66.6666666666667 perc. fin"
## [1] "100 perc. fin"
## [1] "Successfully Finished."
```

``` r
  Result_HLA1_SNV <- CalculatePriorityScores(result = Result_HLA1_SNV, useRNAvaf = FALSE)
  print(head(Result_HLA1_SNV))
```

```
##      HLA           Pos Gene      Evaluated_Mutant_Peptide Mut_EL      Mut_Rank
## [1,] "HLA-A*02:01" "2" "0_DHX15" "ERDYLEAAIRA"            "0.0001110" "39.889"
## [2,] "HLA-A*02:01" "3" "0_DHX15" "RDYLEAAIRAV"            "0.0093530" "6.946" 
## [3,] "HLA-A*02:01" "4" "0_DHX15" "DYLEAAIRAVI"            "0.0003140" "27.742"
## [4,] "HLA-A*02:01" "5" "0_DHX15" "YLEAAIRAVIQ"            "0.0004340" "24.625"
## [5,] "HLA-A*02:01" "6" "0_DHX15" "LEAAIRAVIQI"            "0.0001100" "40.000"
## [6,] "HLA-A*02:01" "7" "0_DHX15" "EAAIRAVIQIH"            "0.0000070" "78.333"
##      Evaluated_Wt_Peptide Wt_EL       Wt_Rank  Chr NM_ID       Change     Ref
## [1,] "ERDYLEAAIRT"        "0.0000310" "58.182" "4" "NM_001358" "c.A1012G" "T"
## [2,] "RDYLEAAIRTV"        "0.0192400" "4.851"  "4" "NM_001358" "c.A1012G" "T"
## [3,] "DYLEAAIRTVI"        "0.0003800" "25.878" "4" "NM_001358" "c.A1012G" "T"
## [4,] "YLEAAIRTVIQ"        "0.0009320" "18.339" "4" "NM_001358" "c.A1012G" "T"
## [5,] "LEAAIRTVIQI"        "0.0001630" "35.077" "4" "NM_001358" "c.A1012G" "T"
## [6,] "EAAIRTVIQIH"        "0.0000080" "76.667" "4" "NM_001358" "c.A1012G" "T"
##      Alt Prob Mutation_Prob. Exon_Start Exon_End   Mutation_Position
## [1,] "C" "0"  "0"            "24529097" "24586177" "4_24556416"     
## [2,] "C" "0"  "0"            "24529097" "24586177" "4_24556416"     
## [3,] "C" "0"  "0"            "24529097" "24586177" "4_24556416"     
## [4,] "C" "0"  "0"            "24529097" "24586177" "4_24556416"     
## [5,] "C" "0"  "0"            "24529097" "24586177" "4_24556416"     
## [6,] "C" "0"  "0"            "24529097" "24586177" "4_24556416"     
##      Total_Depth Tumor_Depth Wt_Peptide               
## [1,] "294"       "143"       "PERDYLEAAIRTVIQIHMCEEEE"
## [2,] "294"       "143"       "PERDYLEAAIRTVIQIHMCEEEE"
## [3,] "294"       "143"       "PERDYLEAAIRTVIQIHMCEEEE"
## [4,] "294"       "143"       "PERDYLEAAIRTVIQIHMCEEEE"
## [5,] "294"       "143"       "PERDYLEAAIRTVIQIHMCEEEE"
## [6,] "294"       "143"       "PERDYLEAAIRTVIQIHMCEEEE"
##      Mutant_Peptide            Total_RNA Tumor_RNA_Ratio Tumor_RNA
## [1,] "PERDYLEAAIRAVIQIHMCEEEE" "1.35204" "NA"            "NA"     
## [2,] "PERDYLEAAIRAVIQIHMCEEEE" "1.35204" "NA"            "NA"     
## [3,] "PERDYLEAAIRAVIQIHMCEEEE" "1.35204" "NA"            "NA"     
## [4,] "PERDYLEAAIRAVIQIHMCEEEE" "1.35204" "NA"            "NA"     
## [5,] "PERDYLEAAIRAVIQIHMCEEEE" "1.35204" "NA"            "NA"     
## [6,] "PERDYLEAAIRAVIQIHMCEEEE" "1.35204" "NA"            "NA"     
##      Tumor_RNA_based_on_DNA MutRatio MutRatio_Min MutRatio_Max
## [1,] "0.657624897959184"    "NA"     "NA"         "NA"        
## [2,] "0.657624897959184"    "NA"     "NA"         "NA"        
## [3,] "0.657624897959184"    "NA"     "NA"         "NA"        
## [4,] "0.657624897959184"    "NA"     "NA"         "NA"        
## [5,] "0.657624897959184"    "NA"     "NA"         "NA"        
## [6,] "0.657624897959184"    "NA"     "NA"         "NA"        
##      P_I                 P_R                     P                      
## [1,] "0.28500874209118"  "3.02670563878825e-83"  "2.25864121245377e-83" 
## [2,] "0.285008765600042" "1.03701641904894e-11"  "7.73860527445574e-12" 
## [3,] "0.285008742435297" "7.20851064695389e-57"  "5.37926087656814e-57" 
## [4,] "0.28500874345477"  "4.22985076083931e-50"  "3.15647320589307e-50" 
## [5,] "0.285008742405142" "1.73754775003521e-83"  "1.29662326806482e-83" 
## [6,] "0.28500874228266"  "1.00203957573045e-166" "7.47759495753465e-167"
```


``` r
  print(Export_Summary_SNV(Input = Result_HLA1_SNV, Mut_Rank_th = 5, Wt_Rank_th = 5))
```

```
##              Num_All_Alteration        Num_Evaluated_Alteration 
##                               3                               3 
## Num_Alteration_Generating_NeoAg                 Num_All_Peptide 
##                               2                             114 
##           Num_Evaluated_Peptide    Num_Peptide_Generating_NeoAg 
##                             113                               3
```


``` r
 Result_HLA1_SNV_vep <- MainSNVClass1(input_vep_format_file = "./../data/sample_vcf.vep.txt",
                                   file_name_in_hla_table = "sample",
                                   hla_file = "./../data/sample_hla_table_c1.txt",
                                   refflat_file  = "./../lib/refFlat.grch37.txt",
                                   refmrna_file = "./../lib/refMrna.grch37.fa",
                                   rnaexp_file = "./../data/sample_rna_exp.txt",
                                   netMHCpan_dir = "./../lib/netMHCpan-4.0/netMHCpan",
                                   multiple_variants = FALSE)
```

```
## [1] "Please Confirm that Reading Start Line is 1"
## [1] "Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\tcDNA_position\tCDS_position\tProtein_position\tAmino_acids\tCodons\tExisting_variation\tExtra"
## [1] "Executing Transformation"
## [1] "HLAtype: A*02:01" "HLAtype: A*32:01" "HLAtype: B*15:17" "HLAtype: B*51:01"
## [5] "HLAtype: C*07:01" "HLAtype: C*15:02"
## [1] "Set chr_column as 1 Chr"
## [1] "Set mutation_start_column as 2 Start"
## [1] "Set mutation_end_column as 3 End"
## [1] "Set mutation_ref_column as 4 Ref"
## [1] "Set mutation_alt_column as 5 Alt"
## [1] "Please Manually Indicate depth_normal_column if you want to use."
## [1] "Set depth_normal_column as NA NA"
## [1] "Please Manually Indicate depth_tumor_column if you want to use."
## [1] "Set depth_tumor_column as NA NA"
## [1] "Set nm_id_column as 10 AAChange.refGene"
## [1] "Please Confirm that Reading Start Line is 1"
## [1] "Chr\tStart\tEnd\tRef\tAlt\tFunc.refGene\tGene.refGene\tGeneDetail.refGene\tExonicFunc.refGene\tAAChange.refGene"
## [1] "Start Analysis: Mutation 1"
## [1] "NM_001282584 - Ref Matched to and NM_ID Attached Information"
## [1] "Peptide Successfully Generated!!"
## [1] "NM_001282582 - Ref Matched to and NM_ID Attached Information"
## [1] "Peptide Successfully Generated!!"
## [1] "NM_032348 - Ref Matched to and NM_ID Attached Information"
## [1] "Peptide Successfully Generated!!"
## [1] "NM_001282585 - Ref Matched to and NM_ID Attached Information"
## [1] "Peptide Successfully Generated!!"
## [1] "NM_001282583 - Ref Matched to and NM_ID Attached Information"
## [1] "Peptide Successfully Generated!!"
## [1] "Start Analysis: Mutation 2"
## [1] "NM_023014 - Ref Matched to and NM_ID Attached Information"
## [1] "Peptide Successfully Generated!!"
## [1] "Start Analysis: Mutation 3"
## [1] "NM_ID NOT Macth, Skip: NM_001405693"
## [1] "NM_ID NOT Macth, Skip: NM_001405684"
## [1] "NM_ID NOT Macth, Skip: NM_001405683"
## [1] "NM_ID NOT Macth, Skip: NM_001405673"
## [1] "NM_ID NOT Macth, Skip: NM_001405695"
## [1] "NM_ID NOT Macth, Skip: NM_001405670"
## [1] "NM_ID NOT Macth, Skip: NM_001405675"
## [1] "NM_ID NOT Macth, Skip: NM_001405686"
## [1] "NM_ID NOT Macth, Skip: NM_001405671"
## [1] "NM_ID NOT Macth, Skip: NM_001405696"
## [1] "NM_ID NOT Macth, Skip: NM_001405687"
## [1] "NM_ID NOT Macth, Skip: NM_001405668"
## [1] "NM_ID NOT Macth, Skip: NM_001405678"
## [1] "NM_ID NOT Macth, Skip: NM_001405680"
## [1] "NM_ID NOT Macth, Skip: NM_001405692"
## [1] "NM_ID NOT Macth, Skip: NM_001405669"
## [1] "NM_ID NOT Macth, Skip: NM_001405674"
## [1] "NM_ID NOT Macth, Skip: NM_001405682"
## [1] "NM_ID NOT Macth, Skip: NM_001405699"
## [1] "NM_ID NOT Macth, Skip: NM_001405700"
## [1] "NM_ID NOT Macth, Skip: NM_001405679"
## [1] "NM_ID NOT Macth, Skip: NM_001405697"
## [1] "NM_ID NOT Macth, Skip: NM_001405698"
## [1] "NM_ID NOT Macth, Skip: NM_001405685"
## [1] "NM_ID NOT Macth, Skip: NM_001405676"
## [1] "NM_ID NOT Macth, Skip: NM_001405681"
## [1] "NM_ID NOT Macth, Skip: NM_001405667"
## [1] "NM_ID NOT Macth, Skip: NM_001405694"
## [1] "NM_ID NOT Macth, Skip: NM_001405672"
## [1] "NM_ID NOT Macth, Skip: NM_001405677"
## [1] "NM_017940 - Ref Matched to and NM_ID Attached Information"
## [1] "End Position Amino Acid is not X, Skip NM_017940"
## [1] -999
## [1] "NM_ID NOT Macth, Skip: NM_001405666"
## [1] "Start Analysis: Mutation 4"
## [1] "NM_ID NOT Macth, Skip: NM_001405693"
## [1] "NM_ID NOT Macth, Skip: NM_001405684"
## [1] "NM_ID NOT Macth, Skip: NM_001405683"
## [1] "NM_ID NOT Macth, Skip: NM_001405673"
## [1] "NM_ID NOT Macth, Skip: NM_001405695"
## [1] "NM_ID NOT Macth, Skip: NM_001405670"
## [1] "NM_ID NOT Macth, Skip: NM_001405675"
## [1] "NM_ID NOT Macth, Skip: NM_001405686"
## [1] "NM_ID NOT Macth, Skip: NM_001405671"
## [1] "NM_ID NOT Macth, Skip: NM_001405696"
## [1] "NM_ID NOT Macth, Skip: NM_001405687"
## [1] "NM_ID NOT Macth, Skip: NM_001405668"
## [1] "NM_ID NOT Macth, Skip: NM_001405678"
## [1] "NM_ID NOT Macth, Skip: NM_001405680"
## [1] "NM_ID NOT Macth, Skip: NM_001405692"
## [1] "NM_ID NOT Macth, Skip: NM_001405669"
## [1] "NM_ID NOT Macth, Skip: NM_001405674"
## [1] "NM_ID NOT Macth, Skip: NM_001405682"
## [1] "NM_ID NOT Macth, Skip: NM_001405699"
## [1] "NM_ID NOT Macth, Skip: NM_001405700"
## [1] "NM_ID NOT Macth, Skip: NM_001405679"
## [1] "NM_ID NOT Macth, Skip: NM_001405697"
## [1] "NM_ID NOT Macth, Skip: NM_001405698"
## [1] "NM_ID NOT Macth, Skip: NM_001405685"
## [1] "NM_ID NOT Macth, Skip: NM_001405676"
## [1] "NM_ID NOT Macth, Skip: NM_001405681"
## [1] "NM_ID NOT Macth, Skip: NM_001405667"
## [1] "NM_ID NOT Macth, Skip: NM_001405694"
## [1] "NM_ID NOT Macth, Skip: NM_001405672"
## [1] "NM_ID NOT Macth, Skip: NM_001405677"
## [1] "NM_017940 - Ref Matched to and NM_ID Attached Information"
## [1] "End Position Amino Acid is not X, Skip NM_017940"
## [1] -999
## [1] "NM_ID NOT Macth, Skip: NM_001405666"
## [1] "Start Analysis: Mutation 5"
## [1] "NM_001164586 - Ref Matched to and NM_ID Attached Information"
## [1] "Peptide Successfully Generated!!"
## [1] "NM_001367841 - Ref Matched to and NM_ID Attached Information"
## [1] "The Mutation is not between Exon Region, Skip NM_001367841"
## [1] "result.ID.SNV1/sample_vcf.vep.txt.annovar_format.txt.ID_SNV.peptide.txt.list.vcf"
## [1] "Did not find ./../lib/netMHCpan-4.0/netMHCpan"
```

``` r
  Result_HLA1_vep_SNV <- CalculatePriorityScores(result = Result_HLA1_SNV_vep, useRNAvaf = FALSE)  
  print(head(Result_HLA1_vep_SNV))
```

```
##      result P_I P_R P
```


``` r
  print(Export_Summary_SNV(Input = Result_HLA1_vep_SNV, Mut_Rank_th = 5, Wt_Rank_th = 5))
```

```
##              Num_All_Alteration        Num_Evaluated_Alteration 
##                               0                               0 
## Num_Alteration_Generating_NeoAg                 Num_All_Peptide 
##                               0                               0 
##           Num_Evaluated_Peptide    Num_Peptide_Generating_NeoAg 
##                               0                               0
```

To calculate the binding affinity of neoantigen candaites, which are generated from from SNVs, to HLA ClassII.


``` r
 Result_HLA2_SNV <- MainSNVClass2(input_annovar_format_file = "./../data/sample_vcf.annovar.txt",
                                   file_name_in_hla_table = "sample",
                                   hla_file = "./../data/sample_hla_table_c2.txt",
                                   refflat_file  = "./../lib/refFlat.grch37.txt",
                                   refmrna_file = "./../lib/refMrna.grch37.fa",
                                   rnaexp_file = "./../data/sample_rna_exp.txt",
                                   netMHCIIpan_dir = "./../lib/netMHCIIpan-4.3/netMHCIIpan",
                                   depth_tumor_column = 12,
                                   depth_normal_column = 14,
                                   SNPs = "./../data/sample_vcf.snps.txt",
                                   multiple_variants = TRUE)
```

```
##  [1] "HLAtype: DPA1*01:03" "HLAtype: DPA1*02:01" "HLAtype: DPB1*02:01"
##  [4] "HLAtype: DPB1*09:01" "HLAtype: DQA1*01:02" "HLAtype: DQA1*05:05"
##  [7] "HLAtype: DQB1*03:01" "HLAtype: DQB1*06:04" "HLAtype: DRB1*11:04"
## [10] "HLAtype: DRB1*13:02"
## [1] "Set chr_column as 1 Chr"
## [1] "Set mutation_start_column as 2 Start"
## [1] "Set mutation_end_column as 3 End"
## [1] "Set mutation_ref_column as 4 Ref"
## [1] "Set mutation_alt_column as 5 Alt"
## [1] "Set depth_normal_column as 14 depth_normal"
## [1] "Set depth_tumor_column as 12 depth_tumor"
## [1] "Set nm_id_column as 10 AAChange.refGene"
## [1] "Please Confirm that Reading Start Line is 1"
## [1] "Chr\tStart\tEnd\tRef\tAlt\tFunc.refGene\tGene.refGene\tGeneDetail.refGene\tExonicFunc.refGene\tAAChange.refGene\tcytoBand\tdepth_tumor\tvariantNum_tumor\tdepth_normal"
## [1] "Please Confirm that Reading Start Line is 1"
## [1] "CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT"
## [1] "Start Analysis: Mutation 1"
## [1] "NM_001358 - Ref Matched to and NM_ID Attached Information"
## [1] "Peptide Successfully Generated!!"
## [1] "Start Analysis: Mutation 2"
## [1] "NM_002048 - Ref Matched to and NM_ID Attached Information"
## [1] "Peptide Successfully Generated!!"
## [1] "Start Analysis: Mutation 3"
## [1] "NM_006205 - Ref Matched to and NM_ID Attached Information"
## [1] "Peptide Successfully Generated!!"
## [1] "result.ID.SNV2/sample_vcf.annovar.txt.ID_SNV.peptide.txt.list.vcf"
## [1] "Executing netMHCpan to result.ID.SNV2"
## [1] "./../lib/netMHCIIpan-4.3/netMHCIIpan -length 15 -f result.ID.SNV2/sample_vcf.annovar.txt.ID_SNV.peptide.fasta -choose -cha DPA10103 -choose -chb DPB10201 > result.ID.SNV2/sample_vcf.annovar.txt.ID_SNV.HLACLASS2.1.peptide.txt"
## [1] "./../lib/netMHCIIpan-4.3/netMHCIIpan -length 15 -f result.ID.SNV2/sample_vcf.annovar.txt.ID_SNV.peptide.fasta -choose -cha DPA10103 -choose -chb DPB10901 > result.ID.SNV2/sample_vcf.annovar.txt.ID_SNV.HLACLASS2.2.peptide.txt"
## [1] "./../lib/netMHCIIpan-4.3/netMHCIIpan -length 15 -f result.ID.SNV2/sample_vcf.annovar.txt.ID_SNV.peptide.fasta -choose -cha DPA10201 -choose -chb DPB10201 > result.ID.SNV2/sample_vcf.annovar.txt.ID_SNV.HLACLASS2.3.peptide.txt"
## [1] "./../lib/netMHCIIpan-4.3/netMHCIIpan -length 15 -f result.ID.SNV2/sample_vcf.annovar.txt.ID_SNV.peptide.fasta -choose -cha DPA10201 -choose -chb DPB10901 > result.ID.SNV2/sample_vcf.annovar.txt.ID_SNV.HLACLASS2.4.peptide.txt"
## [1] "./../lib/netMHCIIpan-4.3/netMHCIIpan -length 15 -f result.ID.SNV2/sample_vcf.annovar.txt.ID_SNV.peptide.fasta -choose -cha DQA10102 -choose -chb DQB10301 > result.ID.SNV2/sample_vcf.annovar.txt.ID_SNV.HLACLASS2.5.peptide.txt"
## [1] "./../lib/netMHCIIpan-4.3/netMHCIIpan -length 15 -f result.ID.SNV2/sample_vcf.annovar.txt.ID_SNV.peptide.fasta -choose -cha DQA10102 -choose -chb DQB10604 > result.ID.SNV2/sample_vcf.annovar.txt.ID_SNV.HLACLASS2.6.peptide.txt"
## [1] "./../lib/netMHCIIpan-4.3/netMHCIIpan -length 15 -f result.ID.SNV2/sample_vcf.annovar.txt.ID_SNV.peptide.fasta -choose -cha DQA10505 -choose -chb DQB10301 > result.ID.SNV2/sample_vcf.annovar.txt.ID_SNV.HLACLASS2.7.peptide.txt"
## [1] "./../lib/netMHCIIpan-4.3/netMHCIIpan -length 15 -f result.ID.SNV2/sample_vcf.annovar.txt.ID_SNV.peptide.fasta -choose -cha DQA10505 -choose -chb DQB10604 > result.ID.SNV2/sample_vcf.annovar.txt.ID_SNV.HLACLASS2.8.peptide.txt"
## [1] "./../lib/netMHCIIpan-4.3/netMHCIIpan -length 15 -f result.ID.SNV2/sample_vcf.annovar.txt.ID_SNV.peptide.fasta -a DRB1_1104 > result.ID.SNV2/sample_vcf.annovar.txt.ID_SNV.HLACLASS2.9.peptide.txt"
## [1] "./../lib/netMHCIIpan-4.3/netMHCIIpan -length 15 -f result.ID.SNV2/sample_vcf.annovar.txt.ID_SNV.peptide.fasta -a DRB1_1302 > result.ID.SNV2/sample_vcf.annovar.txt.ID_SNV.HLACLASS2.10.peptide.txt"
## [1] "./../lib/netMHCIIpan-4.3/netMHCIIpan -length 15 -f result.ID.SNV2/sample_vcf.annovar.txt.ID_SNV.wtpeptide.fasta -choose -cha DPA10103 -choose -chb DPB10201 > result.ID.SNV2/sample_vcf.annovar.txt.ID_SNV.HLACLASS2.1.wtpeptide.txt"
## [1] "./../lib/netMHCIIpan-4.3/netMHCIIpan -length 15 -f result.ID.SNV2/sample_vcf.annovar.txt.ID_SNV.wtpeptide.fasta -choose -cha DPA10103 -choose -chb DPB10901 > result.ID.SNV2/sample_vcf.annovar.txt.ID_SNV.HLACLASS2.2.wtpeptide.txt"
## [1] "./../lib/netMHCIIpan-4.3/netMHCIIpan -length 15 -f result.ID.SNV2/sample_vcf.annovar.txt.ID_SNV.wtpeptide.fasta -choose -cha DPA10201 -choose -chb DPB10201 > result.ID.SNV2/sample_vcf.annovar.txt.ID_SNV.HLACLASS2.3.wtpeptide.txt"
## [1] "./../lib/netMHCIIpan-4.3/netMHCIIpan -length 15 -f result.ID.SNV2/sample_vcf.annovar.txt.ID_SNV.wtpeptide.fasta -choose -cha DPA10201 -choose -chb DPB10901 > result.ID.SNV2/sample_vcf.annovar.txt.ID_SNV.HLACLASS2.4.wtpeptide.txt"
## [1] "./../lib/netMHCIIpan-4.3/netMHCIIpan -length 15 -f result.ID.SNV2/sample_vcf.annovar.txt.ID_SNV.wtpeptide.fasta -choose -cha DQA10102 -choose -chb DQB10301 > result.ID.SNV2/sample_vcf.annovar.txt.ID_SNV.HLACLASS2.5.wtpeptide.txt"
## [1] "./../lib/netMHCIIpan-4.3/netMHCIIpan -length 15 -f result.ID.SNV2/sample_vcf.annovar.txt.ID_SNV.wtpeptide.fasta -choose -cha DQA10102 -choose -chb DQB10604 > result.ID.SNV2/sample_vcf.annovar.txt.ID_SNV.HLACLASS2.6.wtpeptide.txt"
## [1] "./../lib/netMHCIIpan-4.3/netMHCIIpan -length 15 -f result.ID.SNV2/sample_vcf.annovar.txt.ID_SNV.wtpeptide.fasta -choose -cha DQA10505 -choose -chb DQB10301 > result.ID.SNV2/sample_vcf.annovar.txt.ID_SNV.HLACLASS2.7.wtpeptide.txt"
## [1] "./../lib/netMHCIIpan-4.3/netMHCIIpan -length 15 -f result.ID.SNV2/sample_vcf.annovar.txt.ID_SNV.wtpeptide.fasta -choose -cha DQA10505 -choose -chb DQB10604 > result.ID.SNV2/sample_vcf.annovar.txt.ID_SNV.HLACLASS2.8.wtpeptide.txt"
## [1] "./../lib/netMHCIIpan-4.3/netMHCIIpan -length 15 -f result.ID.SNV2/sample_vcf.annovar.txt.ID_SNV.wtpeptide.fasta -a DRB1_1104 > result.ID.SNV2/sample_vcf.annovar.txt.ID_SNV.HLACLASS2.9.wtpeptide.txt"
## [1] "./../lib/netMHCIIpan-4.3/netMHCIIpan -length 15 -f result.ID.SNV2/sample_vcf.annovar.txt.ID_SNV.wtpeptide.fasta -a DRB1_1302 > result.ID.SNV2/sample_vcf.annovar.txt.ID_SNV.HLACLASS2.10.wtpeptide.txt"
## [1] "Merging Results..."
## [1] "/Users/takanorihasegawa/Documents/GitHub/Neoantimon/vignettes/result.ID.SNV2/sample_vcf.annovar.txt.ID_SNV.HLACLASS2.1.peptide.txt"
## [1] "33.3333333333333 perc. fin"
## [1] "66.6666666666667 perc. fin"
## [1] "100 perc. fin"
## [1] "/Users/takanorihasegawa/Documents/GitHub/Neoantimon/vignettes/result.ID.SNV2/sample_vcf.annovar.txt.ID_SNV.HLACLASS2.10.peptide.txt"
## [1] "33.3333333333333 perc. fin"
## [1] "66.6666666666667 perc. fin"
## [1] "100 perc. fin"
## [1] "/Users/takanorihasegawa/Documents/GitHub/Neoantimon/vignettes/result.ID.SNV2/sample_vcf.annovar.txt.ID_SNV.HLACLASS2.2.peptide.txt"
## [1] "33.3333333333333 perc. fin"
## [1] "66.6666666666667 perc. fin"
## [1] "100 perc. fin"
## [1] "/Users/takanorihasegawa/Documents/GitHub/Neoantimon/vignettes/result.ID.SNV2/sample_vcf.annovar.txt.ID_SNV.HLACLASS2.3.peptide.txt"
## [1] "33.3333333333333 perc. fin"
## [1] "66.6666666666667 perc. fin"
## [1] "100 perc. fin"
## [1] "/Users/takanorihasegawa/Documents/GitHub/Neoantimon/vignettes/result.ID.SNV2/sample_vcf.annovar.txt.ID_SNV.HLACLASS2.4.peptide.txt"
## [1] "33.3333333333333 perc. fin"
## [1] "66.6666666666667 perc. fin"
## [1] "100 perc. fin"
## [1] "/Users/takanorihasegawa/Documents/GitHub/Neoantimon/vignettes/result.ID.SNV2/sample_vcf.annovar.txt.ID_SNV.HLACLASS2.5.peptide.txt"
## [1] "33.3333333333333 perc. fin"
## [1] "66.6666666666667 perc. fin"
## [1] "100 perc. fin"
## [1] "/Users/takanorihasegawa/Documents/GitHub/Neoantimon/vignettes/result.ID.SNV2/sample_vcf.annovar.txt.ID_SNV.HLACLASS2.6.peptide.txt"
## [1] "33.3333333333333 perc. fin"
## [1] "66.6666666666667 perc. fin"
## [1] "100 perc. fin"
## [1] "/Users/takanorihasegawa/Documents/GitHub/Neoantimon/vignettes/result.ID.SNV2/sample_vcf.annovar.txt.ID_SNV.HLACLASS2.7.peptide.txt"
## [1] "33.3333333333333 perc. fin"
## [1] "66.6666666666667 perc. fin"
## [1] "100 perc. fin"
## [1] "/Users/takanorihasegawa/Documents/GitHub/Neoantimon/vignettes/result.ID.SNV2/sample_vcf.annovar.txt.ID_SNV.HLACLASS2.8.peptide.txt"
## [1] "33.3333333333333 perc. fin"
## [1] "66.6666666666667 perc. fin"
## [1] "100 perc. fin"
## [1] "/Users/takanorihasegawa/Documents/GitHub/Neoantimon/vignettes/result.ID.SNV2/sample_vcf.annovar.txt.ID_SNV.HLACLASS2.9.peptide.txt"
## [1] "33.3333333333333 perc. fin"
## [1] "66.6666666666667 perc. fin"
## [1] "100 perc. fin"
## [1] "Successfully Finished."
```

``` r
  Result_HLA2_SNV <- CalculatePriorityScores(result = Result_HLA2_SNV, useRNAvaf = FALSE) 
  print(head(Result_HLA2_SNV))
```

```
##      HLA         Pos Gene      Evaluated_Mutant_Peptide Mut_EL     Mut_Rank
## [1,] "DRB1_0101" "2" "0_DHX15" "TPEPERDYLEAAIRA"        "0.000549" "55.59" 
## [2,] "DRB1_0101" "3" "0_DHX15" "PEPERDYLEAAIRAV"        "0.007817" "26.45" 
## [3,] "DRB1_0101" "4" "0_DHX15" "EPERDYLEAAIRAVI"        "0.036515" "14.73" 
## [4,] "DRB1_0101" "5" "0_DHX15" "PERDYLEAAIRAVIQ"        "0.113296" "8.52"  
## [5,] "DRB1_0101" "6" "0_DHX15" "ERDYLEAAIRAVIQI"        "0.093106" "9.46"  
## [6,] "DRB1_0101" "7" "0_DHX15" "RDYLEAAIRAVIQIH"        "0.012420" "22.51" 
##      Evaluated_Wt_Peptide Wt_EL      Wt_Rank Chr NM_ID       Change     Ref Alt
## [1,] "TPEPERDYLEAAIRT"    "0.000295" "63.70" "4" "NM_001358" "c.A1012G" "T" "C"
## [2,] "PEPERDYLEAAIRTV"    "0.004001" "32.76" "4" "NM_001358" "c.A1012G" "T" "C"
## [3,] "EPERDYLEAAIRTVI"    "0.027292" "16.63" "4" "NM_001358" "c.A1012G" "T" "C"
## [4,] "PERDYLEAAIRTVIQ"    "0.082495" "10.06" "4" "NM_001358" "c.A1012G" "T" "C"
## [5,] "ERDYLEAAIRTVIQI"    "0.077206" "10.41" "4" "NM_001358" "c.A1012G" "T" "C"
## [6,] "RDYLEAAIRTVIQIH"    "0.011714" "22.96" "4" "NM_001358" "c.A1012G" "T" "C"
##      Prob Mutation_Prob. Exon_Start Exon_End   Mutation_Position Total_Depth
## [1,] "0"  "0"            "24529097" "24586177" "4_24556416"      "294"      
## [2,] "0"  "0"            "24529097" "24586177" "4_24556416"      "294"      
## [3,] "0"  "0"            "24529097" "24586177" "4_24556416"      "294"      
## [4,] "0"  "0"            "24529097" "24586177" "4_24556416"      "294"      
## [5,] "0"  "0"            "24529097" "24586177" "4_24556416"      "294"      
## [6,] "0"  "0"            "24529097" "24586177" "4_24556416"      "294"      
##      Tumor_Depth Wt_Peptide                       
## [1,] "143"       "YTPEPERDYLEAAIRTVIQIHMCEEEEGDLL"
## [2,] "143"       "YTPEPERDYLEAAIRTVIQIHMCEEEEGDLL"
## [3,] "143"       "YTPEPERDYLEAAIRTVIQIHMCEEEEGDLL"
## [4,] "143"       "YTPEPERDYLEAAIRTVIQIHMCEEEEGDLL"
## [5,] "143"       "YTPEPERDYLEAAIRTVIQIHMCEEEEGDLL"
## [6,] "143"       "YTPEPERDYLEAAIRTVIQIHMCEEEEGDLL"
##      Mutant_Peptide                    Total_RNA Tumor_RNA_Ratio Tumor_RNA
## [1,] "YTPEPERDYLEAAIRAVIQIHMCEEEEGDLL" "1.35204" "NA"            "NA"     
## [2,] "YTPEPERDYLEAAIRAVIQIHMCEEEEGDLL" "1.35204" "NA"            "NA"     
## [3,] "YTPEPERDYLEAAIRAVIQIHMCEEEEGDLL" "1.35204" "NA"            "NA"     
## [4,] "YTPEPERDYLEAAIRAVIQIHMCEEEEGDLL" "1.35204" "NA"            "NA"     
## [5,] "YTPEPERDYLEAAIRAVIQIHMCEEEEGDLL" "1.35204" "NA"            "NA"     
## [6,] "YTPEPERDYLEAAIRAVIQIHMCEEEEGDLL" "1.35204" "NA"            "NA"     
##      Tumor_RNA_based_on_DNA MutRatio MutRatio_Min MutRatio_Max
## [1,] "0.657624897959184"    "NA"     "NA"         "NA"        
## [2,] "0.657624897959184"    "NA"     "NA"         "NA"        
## [3,] "0.657624897959184"    "NA"     "NA"         "NA"        
## [4,] "0.657624897959184"    "NA"     "NA"         "NA"        
## [5,] "0.657624897959184"    "NA"     "NA"         "NA"        
## [6,] "0.657624897959184"    "NA"     "NA"         "NA"        
##      P_I                 P_R                     P                     
## [1,] "0.28500874167929"  "2.43602580240382e-117" "1.8178537752065e-117"
## [2,] "0.285008733251095" "4.60667843033526e-54"  "3.43767613934292e-54"
## [3,] "0.285008720402792" "1.2973123034524e-28"   "9.68103074329359e-29"
## [4,] "0.285008669168784" "3.96176033890512e-15"  "2.95641408290313e-15"
## [5,] "0.285008704455587" "3.60333080928503e-17"  "2.68894053113598e-17"
## [6,] "0.285008740581003" "1.65572875348136e-45"  "1.23556686561528e-45"
```


``` r
  print(Export_Summary_SNV(Input = Result_HLA2_SNV, Mut_Rank_th = 10, Wt_Rank_th = 10))
```

```
##              Num_All_Alteration        Num_Evaluated_Alteration 
##                               3                               3 
## Num_Alteration_Generating_NeoAg                 Num_All_Peptide 
##                               1                              45 
##           Num_Evaluated_Peptide    Num_Peptide_Generating_NeoAg 
##                              45                               3
```

To calculate the binding affinity of neoantigen candaites, which are generated from from indels, to HLA ClassI. When using MHCflurry, [[1]] and [[2]] include the results of NetMHCpan and MHCflurry, respectively.


``` r
Result_HLA1_INDEL <- MainINDELClass1(input_annovar_format_file = "./../data/sample_vcf.annovar.txt",
                                       file_name_in_hla_table = "sample",
                                       hla_file = "./../data/sample_hla_table_c1.txt",
                                       refflat_file  = "./../lib/refFlat.grch37.txt",
                                       refmrna_file = "./../lib/refMrna.grch37.fa",
                                       rnaexp_file = "./../data/sample_rna_exp.txt",
                                       netMHCpan_dir = "./../lib/netMHCpan-4.1/netMHCpan",
                                       depth_tumor_column = 12,
                                       depth_normal_column = 14,
                                       SNPs = "./../data/sample_vcf.snps.txt",
                                       multiple_variants = TRUE)
```

```
## [1] "HLAtype: A*02:01" "HLAtype: A*32:01" "HLAtype: B*15:17" "HLAtype: B*51:01"
## [5] "HLAtype: C*07:01" "HLAtype: C*15:02"
## [1] "Set chr_column as 1 Chr"
## [1] "Set mutation_start_column as 2 Start"
## [1] "Set mutation_end_column as 3 End"
## [1] "Set mutation_ref_column as 4 Ref"
## [1] "Set mutation_alt_column as 5 Alt"
## [1] "Set depth_normal_column as 14 depth_normal"
## [1] "Set depth_tumor_column as 12 depth_tumor"
## [1] "Set nm_id_column as 10 AAChange.refGene"
## [1] "Please Confirm that Reading Start Line is 1"
## [1] "Chr\tStart\tEnd\tRef\tAlt\tFunc.refGene\tGene.refGene\tGeneDetail.refGene\tExonicFunc.refGene\tAAChange.refGene\tcytoBand\tdepth_tumor\tvariantNum_tumor\tdepth_normal"
## [1] "Please Confirm that Reading Start Line is 1"
## [1] "CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT"
## [1] "Start Analysis: Mutation 1"
## [1] "Peptide Successfully Generated!!"
## [1] "Start Analysis: Mutation 2"
## [1] "Peptide Successfully Generated!!"
## [1] "Start Analysis: Mutation 3"
## [1] "Peptide Successfully Generated!!"
## [1] "Peptide Successfully Generated!!"
## [1] "Peptide Successfully Generated!!"
## [1] "Peptide Successfully Generated!!"
## [1] "result.ID.INDEL1/sample_vcf.annovar.txt.ID_INDEL.peptide.txt.list.vcf"
## [1] "Executing netMHCpan to result.ID.INDEL1"
## [1] "./../lib/netMHCpan-4.1/netMHCpan -BA  -l 8,9,10,11 -f result.ID.INDEL1/sample_vcf.annovar.txt.ID_INDEL.peptide.fasta -a HLA-A02:01 > result.ID.INDEL1/sample_vcf.annovar.txt.ID_INDEL.HLACLASS1.1.peptide.txt"
## [1] "./../lib/netMHCpan-4.1/netMHCpan -BA  -l 8,9,10,11 -f result.ID.INDEL1/sample_vcf.annovar.txt.ID_INDEL.peptide.fasta -a HLA-A32:01 > result.ID.INDEL1/sample_vcf.annovar.txt.ID_INDEL.HLACLASS1.2.peptide.txt"
## [1] "./../lib/netMHCpan-4.1/netMHCpan -BA  -l 8,9,10,11 -f result.ID.INDEL1/sample_vcf.annovar.txt.ID_INDEL.peptide.fasta -a HLA-B15:17 > result.ID.INDEL1/sample_vcf.annovar.txt.ID_INDEL.HLACLASS1.3.peptide.txt"
## [1] "./../lib/netMHCpan-4.1/netMHCpan -BA  -l 8,9,10,11 -f result.ID.INDEL1/sample_vcf.annovar.txt.ID_INDEL.peptide.fasta -a HLA-B51:01 > result.ID.INDEL1/sample_vcf.annovar.txt.ID_INDEL.HLACLASS1.4.peptide.txt"
## [1] "./../lib/netMHCpan-4.1/netMHCpan -BA  -l 8,9,10,11 -f result.ID.INDEL1/sample_vcf.annovar.txt.ID_INDEL.peptide.fasta -a HLA-C07:01 > result.ID.INDEL1/sample_vcf.annovar.txt.ID_INDEL.HLACLASS1.5.peptide.txt"
## [1] "./../lib/netMHCpan-4.1/netMHCpan -BA  -l 8,9,10,11 -f result.ID.INDEL1/sample_vcf.annovar.txt.ID_INDEL.peptide.fasta -a HLA-C15:02 > result.ID.INDEL1/sample_vcf.annovar.txt.ID_INDEL.HLACLASS1.6.peptide.txt"
## [1] "Merging Results..."
## [1] "/Users/takanorihasegawa/Documents/GitHub/Neoantimon/vignettes/result.ID.INDEL1/sample_vcf.annovar.txt.ID_INDEL.HLACLASS1.1.peptide.txt"
## [1] "33.3333333333333 perc. fin"
## [1] "66.6666666666667 perc. fin"
## [1] "100 perc. fin"
## [1] "/Users/takanorihasegawa/Documents/GitHub/Neoantimon/vignettes/result.ID.INDEL1/sample_vcf.annovar.txt.ID_INDEL.HLACLASS1.2.peptide.txt"
## [1] "33.3333333333333 perc. fin"
## [1] "66.6666666666667 perc. fin"
## [1] "100 perc. fin"
## [1] "/Users/takanorihasegawa/Documents/GitHub/Neoantimon/vignettes/result.ID.INDEL1/sample_vcf.annovar.txt.ID_INDEL.HLACLASS1.3.peptide.txt"
## [1] "33.3333333333333 perc. fin"
## [1] "66.6666666666667 perc. fin"
## [1] "100 perc. fin"
## [1] "/Users/takanorihasegawa/Documents/GitHub/Neoantimon/vignettes/result.ID.INDEL1/sample_vcf.annovar.txt.ID_INDEL.HLACLASS1.4.peptide.txt"
## [1] "33.3333333333333 perc. fin"
## [1] "66.6666666666667 perc. fin"
## [1] "100 perc. fin"
## [1] "/Users/takanorihasegawa/Documents/GitHub/Neoantimon/vignettes/result.ID.INDEL1/sample_vcf.annovar.txt.ID_INDEL.HLACLASS1.5.peptide.txt"
## [1] "33.3333333333333 perc. fin"
## [1] "66.6666666666667 perc. fin"
## [1] "100 perc. fin"
## [1] "/Users/takanorihasegawa/Documents/GitHub/Neoantimon/vignettes/result.ID.INDEL1/sample_vcf.annovar.txt.ID_INDEL.HLACLASS1.6.peptide.txt"
## [1] "33.3333333333333 perc. fin"
## [1] "66.6666666666667 perc. fin"
## [1] "100 perc. fin"
## [1] "Successfully Finished."
```

``` r
  Result_HLA1_INDEL_1 <- CalculatePriorityScores(result = Result_HLA1_INDEL, useRNAvaf = FALSE)
  print(head(Result_HLA1_INDEL_1))
```

```
##      HLA           Pos Gene        Evaluated_Mutant_Peptide_Core
## [1,] "HLA-A*02:01" "1" "0_UGT2B28" "YHGIPMVGIPF"                
## [2,] "HLA-A*02:01" "2" "0_UGT2B28" "HGIPMVGIPFV"                
## [3,] "HLA-A*02:01" "3" "0_UGT2B28" "GIPMVGIPFVL"                
## [4,] "HLA-A*02:01" "4" "0_UGT2B28" "IPMVGIPFVL"                 
## [5,] "HLA-A*02:01" "5" "0_UGT2B28" "PMVGIPFVLGS"                
## [6,] "HLA-A*02:01" "6" "0_UGT2B28" "MVGIPFVLGST"                
##      Evaluated_Mutant_Peptide Mut_EL      Mut_Rank Chr NM_ID      
## [1,] "YHGIPMVGIPF"            "0.0000510" "50.667" "4" "NM_053039"
## [2,] "HGIPMVGIPFV"            "0.0034220" "10.843" "4" "NM_053039"
## [3,] "GIPMVGIPFVL"            "0.0346280" "3.599"  "4" "NM_053039"
## [4,] "IPMVGIPFVLG"            "0.0006810" "20.759" "4" "NM_053039"
## [5,] "PMVGIPFVLGS"            "0.0000390" "54.667" "4" "NM_053039"
## [6,] "MVGIPFVLGST"            "0.0001620" "35.154" "4" "NM_053039"
##      Change           Ref Alt Prob Mutation_Prob. Exon_Start Exon_End  
## [1,] "Out_c.1186dupT" "-" "T" "0"  "0"            "70146192" "70160768"
## [2,] "Out_c.1186dupT" "-" "T" "0"  "0"            "70146192" "70160768"
## [3,] "Out_c.1186dupT" "-" "T" "0"  "0"            "70146192" "70160768"
## [4,] "Out_c.1186dupT" "-" "T" "0"  "0"            "70146192" "70160768"
## [5,] "Out_c.1186dupT" "-" "T" "0"  "0"            "70146192" "70160768"
## [6,] "Out_c.1186dupT" "-" "T" "0"  "0"            "70146192" "70160768"
##      Mutation_Position Total_Depth Tumor_Depth Wt_Peptide                   
## [1,] "4_70156404"      "84"        "43"        "YHGIPMVGIPLFWDQPDNIAHMKAKGA"
## [2,] "4_70156404"      "84"        "43"        "YHGIPMVGIPLFWDQPDNIAHMKAKGA"
## [3,] "4_70156404"      "84"        "43"        "YHGIPMVGIPLFWDQPDNIAHMKAKGA"
## [4,] "4_70156404"      "84"        "43"        "YHGIPMVGIPLFWDQPDNIAHMKAKGA"
## [5,] "4_70156404"      "84"        "43"        "YHGIPMVGIPLFWDQPDNIAHMKAKGA"
## [6,] "4_70156404"      "84"        "43"        "YHGIPMVGIPLFWDQPDNIAHMKAKGA"
##      Mutant_Peptide      Total_RNA Tumor_RNA_Ratio Tumor_RNA
## [1,] "YHGIPMVGIPFVLGSTX" "0"       "NA"            "NA"     
## [2,] "YHGIPMVGIPFVLGSTX" "0"       "NA"            "NA"     
## [3,] "YHGIPMVGIPFVLGSTX" "0"       "NA"            "NA"     
## [4,] "YHGIPMVGIPFVLGSTX" "0"       "NA"            "NA"     
## [5,] "YHGIPMVGIPFVLGSTX" "0"       "NA"            "NA"     
## [6,] "YHGIPMVGIPFVLGSTX" "0"       "NA"            "NA"     
##      Tumor_RNA_based_on_DNA MutRatio MutRatio_Min MutRatio_Max P_I P_R P  
## [1,] "0"                    "NA"     "NA"         "NA"         "0" "0" "0"
## [2,] "0"                    "NA"     "NA"         "NA"         "0" "0" "0"
## [3,] "0"                    "NA"     "NA"         "NA"         "0" "0" "0"
## [4,] "0"                    "NA"     "NA"         "NA"         "0" "0" "0"
## [5,] "0"                    "NA"     "NA"         "NA"         "0" "0" "0"
## [6,] "0"                    "NA"     "NA"         "NA"         "0" "0" "0"
```


``` r
  print(Export_Summary_IndelSV(Input = Result_HLA1_INDEL_1, Mut_Rank_th = 5))
```

```
##              Num_All_Alteration        Num_Evaluated_Alteration 
##                               3                               3 
## Num_Alteration_Generating_NeoAg                 Num_All_Peptide 
##                               3                             184 
##           Num_Evaluated_Peptide    Num_Peptide_Generating_NeoAg 
##                             184                              53
```


``` r
  print(Export_Summary_IndelSV_perFragments(Input = Result_HLA1_INDEL_1, Mut_Rank_th = 5))
```

```
##                                                     YHGIPMVGIPFVLGSTX-0_UGT2B28
## Num_Peptide_Per_Pep                                                      24.000
## Num_Cond_Peptide_Per_Pep                                                 24.000
## Num_Rest_Peptide_Per_Pep                                                 10.000
## Num_Rest_Peptide_Per_Pep / Num_Cond_Peptide_Per_Pep                       0.417
## -logP                                                                     0.391
##                                                     RDELLAAIKKFHIKVATQX-1_COL12A1
## Num_Peptide_Per_Pep                                                        32.000
## Num_Cond_Peptide_Per_Pep                                                   32.000
## Num_Rest_Peptide_Per_Pep                                                    7.000
## Num_Rest_Peptide_Per_Pep / Num_Cond_Peptide_Per_Pep                         0.219
## -logP                                                                       0.437
##                                                     ARDFLPSLKNRFWKPSILPIFMYKHCSVQFSVRHGDVQTKVHX-2_SLCO1C1
## Num_Peptide_Per_Pep                                                                               128.000
## Num_Cond_Peptide_Per_Pep                                                                          128.000
## Num_Rest_Peptide_Per_Pep                                                                           36.000
## Num_Rest_Peptide_Per_Pep / Num_Cond_Peptide_Per_Pep                                                 0.281
## -logP                                                                                               0.989
```

To calculate the binding affinity of neoantigen candaites, which are generated from from indels, to HLA ClassII.


``` r
  Result_HLA2_INDEL <- MainINDELClass2(input_annovar_format_file = "./../data/sample_vcf.annovar.txt",
                                       file_name_in_hla_table = "sample",
                                       hla_file = "./../data/sample_hla_table_c2.txt",
                                       refflat_file  = "./../lib/refFlat.grch37.txt",
                                       refmrna_file = "./../lib/refMrna.grch37.fa",
                                       rnaexp_file = "./../data/sample_rna_exp.txt",
                                       netMHCIIpan_dir = "./../lib/netMHCIIpan-4.3/netMHCIIpan",
                                       depth_tumor_column = 12,
                                       depth_normal_column = 14,
                                       SNPs = "./../data/sample_vcf.snps.txt",
                                       multiple_variants = TRUE)
```

```
##  [1] "HLAtype: DPA1*01:03" "HLAtype: DPA1*02:01" "HLAtype: DPB1*02:01"
##  [4] "HLAtype: DPB1*09:01" "HLAtype: DQA1*01:02" "HLAtype: DQA1*05:05"
##  [7] "HLAtype: DQB1*03:01" "HLAtype: DQB1*06:04" "HLAtype: DRB1*11:04"
## [10] "HLAtype: DRB1*13:02"
## [1] "Set chr_column as 1 Chr"
## [1] "Set mutation_start_column as 2 Start"
## [1] "Set mutation_end_column as 3 End"
## [1] "Set mutation_ref_column as 4 Ref"
## [1] "Set mutation_alt_column as 5 Alt"
## [1] "Set depth_normal_column as 14 depth_normal"
## [1] "Set depth_tumor_column as 12 depth_tumor"
## [1] "Set nm_id_column as 10 AAChange.refGene"
## [1] "Please Confirm that Reading Start Line is 1"
## [1] "Chr\tStart\tEnd\tRef\tAlt\tFunc.refGene\tGene.refGene\tGeneDetail.refGene\tExonicFunc.refGene\tAAChange.refGene\tcytoBand\tdepth_tumor\tvariantNum_tumor\tdepth_normal"
## [1] "Please Confirm that Reading Start Line is 1"
## [1] "CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT"
## [1] "Start Analysis: Mutation 1"
## [1] "Peptide Successfully Generated!!"
## [1] "Start Analysis: Mutation 2"
## [1] "Peptide Successfully Generated!!"
## [1] "Start Analysis: Mutation 3"
## [1] "Peptide Successfully Generated!!"
## [1] "Peptide Successfully Generated!!"
## [1] "Peptide Successfully Generated!!"
## [1] "Peptide Successfully Generated!!"
## [1] "result.ID.INDEL2/sample_vcf.annovar.txt.ID_INDEL.peptide.txt.list.vcf"
## [1] "Executing netMHCpan to result.ID.INDEL2"
## [1] "./../lib/netMHCIIpan-4.3/netMHCIIpan -length 15 -f result.ID.INDEL2/sample_vcf.annovar.txt.ID_INDEL.peptide.fasta -choose -cha DPA10103 -choose -chb DPB10201 > result.ID.INDEL2/sample_vcf.annovar.txt.ID_INDEL.HLACLASS2.1.peptide.txt"
## [1] "./../lib/netMHCIIpan-4.3/netMHCIIpan -length 15 -f result.ID.INDEL2/sample_vcf.annovar.txt.ID_INDEL.peptide.fasta -choose -cha DPA10103 -choose -chb DPB10901 > result.ID.INDEL2/sample_vcf.annovar.txt.ID_INDEL.HLACLASS2.2.peptide.txt"
## [1] "./../lib/netMHCIIpan-4.3/netMHCIIpan -length 15 -f result.ID.INDEL2/sample_vcf.annovar.txt.ID_INDEL.peptide.fasta -choose -cha DPA10201 -choose -chb DPB10201 > result.ID.INDEL2/sample_vcf.annovar.txt.ID_INDEL.HLACLASS2.3.peptide.txt"
## [1] "./../lib/netMHCIIpan-4.3/netMHCIIpan -length 15 -f result.ID.INDEL2/sample_vcf.annovar.txt.ID_INDEL.peptide.fasta -choose -cha DPA10201 -choose -chb DPB10901 > result.ID.INDEL2/sample_vcf.annovar.txt.ID_INDEL.HLACLASS2.4.peptide.txt"
## [1] "./../lib/netMHCIIpan-4.3/netMHCIIpan -length 15 -f result.ID.INDEL2/sample_vcf.annovar.txt.ID_INDEL.peptide.fasta -choose -cha DQA10102 -choose -chb DQB10301 > result.ID.INDEL2/sample_vcf.annovar.txt.ID_INDEL.HLACLASS2.5.peptide.txt"
## [1] "./../lib/netMHCIIpan-4.3/netMHCIIpan -length 15 -f result.ID.INDEL2/sample_vcf.annovar.txt.ID_INDEL.peptide.fasta -choose -cha DQA10102 -choose -chb DQB10604 > result.ID.INDEL2/sample_vcf.annovar.txt.ID_INDEL.HLACLASS2.6.peptide.txt"
## [1] "./../lib/netMHCIIpan-4.3/netMHCIIpan -length 15 -f result.ID.INDEL2/sample_vcf.annovar.txt.ID_INDEL.peptide.fasta -choose -cha DQA10505 -choose -chb DQB10301 > result.ID.INDEL2/sample_vcf.annovar.txt.ID_INDEL.HLACLASS2.7.peptide.txt"
## [1] "./../lib/netMHCIIpan-4.3/netMHCIIpan -length 15 -f result.ID.INDEL2/sample_vcf.annovar.txt.ID_INDEL.peptide.fasta -choose -cha DQA10505 -choose -chb DQB10604 > result.ID.INDEL2/sample_vcf.annovar.txt.ID_INDEL.HLACLASS2.8.peptide.txt"
## [1] "./../lib/netMHCIIpan-4.3/netMHCIIpan -length 15 -f result.ID.INDEL2/sample_vcf.annovar.txt.ID_INDEL.peptide.fasta -a DRB1_1104 > result.ID.INDEL2/sample_vcf.annovar.txt.ID_INDEL.HLACLASS2.9.peptide.txt"
## [1] "./../lib/netMHCIIpan-4.3/netMHCIIpan -length 15 -f result.ID.INDEL2/sample_vcf.annovar.txt.ID_INDEL.peptide.fasta -a DRB1_1302 > result.ID.INDEL2/sample_vcf.annovar.txt.ID_INDEL.HLACLASS2.10.peptide.txt"
## [1] "Merging Results..."
## [1] "/Users/takanorihasegawa/Documents/GitHub/Neoantimon/vignettes/result.ID.INDEL2/sample_vcf.annovar.txt.ID_INDEL.HLACLASS2.1.peptide.txt"
## [1] "33.3333333333333 perc. fin"
## [1] "66.6666666666667 perc. fin"
## [1] "100 perc. fin"
## [1] "/Users/takanorihasegawa/Documents/GitHub/Neoantimon/vignettes/result.ID.INDEL2/sample_vcf.annovar.txt.ID_INDEL.HLACLASS2.10.peptide.txt"
## [1] "33.3333333333333 perc. fin"
## [1] "66.6666666666667 perc. fin"
## [1] "100 perc. fin"
## [1] "/Users/takanorihasegawa/Documents/GitHub/Neoantimon/vignettes/result.ID.INDEL2/sample_vcf.annovar.txt.ID_INDEL.HLACLASS2.2.peptide.txt"
## [1] "33.3333333333333 perc. fin"
## [1] "66.6666666666667 perc. fin"
## [1] "100 perc. fin"
## [1] "/Users/takanorihasegawa/Documents/GitHub/Neoantimon/vignettes/result.ID.INDEL2/sample_vcf.annovar.txt.ID_INDEL.HLACLASS2.3.peptide.txt"
## [1] "33.3333333333333 perc. fin"
## [1] "66.6666666666667 perc. fin"
## [1] "100 perc. fin"
## [1] "/Users/takanorihasegawa/Documents/GitHub/Neoantimon/vignettes/result.ID.INDEL2/sample_vcf.annovar.txt.ID_INDEL.HLACLASS2.4.peptide.txt"
## [1] "33.3333333333333 perc. fin"
## [1] "66.6666666666667 perc. fin"
## [1] "100 perc. fin"
## [1] "/Users/takanorihasegawa/Documents/GitHub/Neoantimon/vignettes/result.ID.INDEL2/sample_vcf.annovar.txt.ID_INDEL.HLACLASS2.5.peptide.txt"
## [1] "33.3333333333333 perc. fin"
## [1] "66.6666666666667 perc. fin"
## [1] "100 perc. fin"
## [1] "/Users/takanorihasegawa/Documents/GitHub/Neoantimon/vignettes/result.ID.INDEL2/sample_vcf.annovar.txt.ID_INDEL.HLACLASS2.6.peptide.txt"
## [1] "33.3333333333333 perc. fin"
## [1] "66.6666666666667 perc. fin"
## [1] "100 perc. fin"
## [1] "/Users/takanorihasegawa/Documents/GitHub/Neoantimon/vignettes/result.ID.INDEL2/sample_vcf.annovar.txt.ID_INDEL.HLACLASS2.7.peptide.txt"
## [1] "33.3333333333333 perc. fin"
## [1] "66.6666666666667 perc. fin"
## [1] "100 perc. fin"
## [1] "/Users/takanorihasegawa/Documents/GitHub/Neoantimon/vignettes/result.ID.INDEL2/sample_vcf.annovar.txt.ID_INDEL.HLACLASS2.8.peptide.txt"
## [1] "33.3333333333333 perc. fin"
## [1] "66.6666666666667 perc. fin"
## [1] "100 perc. fin"
## [1] "/Users/takanorihasegawa/Documents/GitHub/Neoantimon/vignettes/result.ID.INDEL2/sample_vcf.annovar.txt.ID_INDEL.HLACLASS2.9.peptide.txt"
## [1] "33.3333333333333 perc. fin"
## [1] "66.6666666666667 perc. fin"
## [1] "100 perc. fin"
## [1] "Successfully Finished."
```

``` r
  Result_HLA2_INDEL <- CalculatePriorityScores(result = Result_HLA2_INDEL, useRNAvaf = FALSE)  
  print(head(Result_HLA2_INDEL))
```

```
##      HLA         Pos Gene        Evaluated_Mutant_Peptide_Core
## [1,] "DRB1_0101" "1" "0_UGT2B28" "YHGIPMVGI"                  
## [2,] "DRB1_0101" "2" "0_UGT2B28" "YHGIPMVGI"                  
## [3,] "DRB1_0101" "3" "0_UGT2B28" "YHGIPMVGI"                  
## [4,] "DRB1_0101" "4" "0_UGT2B28" "GIPMVGIPF"                  
## [5,] "DRB1_0101" "5" "0_UGT2B28" "MVGIPFVLG"                  
## [6,] "DRB1_0101" "6" "0_UGT2B28" "MVGIPFVLG"                  
##      Evaluated_Mutant_Peptide Mut_EL     Mut_Rank Chr NM_ID      
## [1,] "YEAIYHGIPMVGIPF"        "0.068310" "11.01"  "4" "NM_053039"
## [2,] "EAIYHGIPMVGIPFV"        "0.064542" "11.34"  "4" "NM_053039"
## [3,] "AIYHGIPMVGIPFVL"        "0.005254" "30.11"  "4" "NM_053039"
## [4,] "IYHGIPMVGIPFVLG"        "0.000447" "58.35"  "4" "NM_053039"
## [5,] "YHGIPMVGIPFVLGS"        "0.000776" "51.47"  "4" "NM_053039"
## [6,] "HGIPMVGIPFVLGST"        "0.001719" "41.82"  "4" "NM_053039"
##      Change           Ref Alt Prob Mutation_Prob. Exon_Start Exon_End  
## [1,] "Out_c.1186dupT" "-" "T" "0"  "0"            "70146192" "70160768"
## [2,] "Out_c.1186dupT" "-" "T" "0"  "0"            "70146192" "70160768"
## [3,] "Out_c.1186dupT" "-" "T" "0"  "0"            "70146192" "70160768"
## [4,] "Out_c.1186dupT" "-" "T" "0"  "0"            "70146192" "70160768"
## [5,] "Out_c.1186dupT" "-" "T" "0"  "0"            "70146192" "70160768"
## [6,] "Out_c.1186dupT" "-" "T" "0"  "0"            "70146192" "70160768"
##      Mutation_Position Total_Depth Tumor_Depth
## [1,] "4_70156404"      "84"        "43"       
## [2,] "4_70156404"      "84"        "43"       
## [3,] "4_70156404"      "84"        "43"       
## [4,] "4_70156404"      "84"        "43"       
## [5,] "4_70156404"      "84"        "43"       
## [6,] "4_70156404"      "84"        "43"       
##      Wt_Peptide                        Mutant_Peptide          Total_RNA
## [1,] "YEAIYHGIPMVGIPLFWDQPDNIAHMKAKGA" "YEAIYHGIPMVGIPFVLGSTX" "0"      
## [2,] "YEAIYHGIPMVGIPLFWDQPDNIAHMKAKGA" "YEAIYHGIPMVGIPFVLGSTX" "0"      
## [3,] "YEAIYHGIPMVGIPLFWDQPDNIAHMKAKGA" "YEAIYHGIPMVGIPFVLGSTX" "0"      
## [4,] "YEAIYHGIPMVGIPLFWDQPDNIAHMKAKGA" "YEAIYHGIPMVGIPFVLGSTX" "0"      
## [5,] "YEAIYHGIPMVGIPLFWDQPDNIAHMKAKGA" "YEAIYHGIPMVGIPFVLGSTX" "0"      
## [6,] "YEAIYHGIPMVGIPLFWDQPDNIAHMKAKGA" "YEAIYHGIPMVGIPFVLGSTX" "0"      
##      Tumor_RNA_Ratio Tumor_RNA Tumor_RNA_based_on_DNA MutRatio MutRatio_Min
## [1,] "NA"            "NA"      "0"                    "NA"     "NA"        
## [2,] "NA"            "NA"      "0"                    "NA"     "NA"        
## [3,] "NA"            "NA"      "0"                    "NA"     "NA"        
## [4,] "NA"            "NA"      "0"                    "NA"     "NA"        
## [5,] "NA"            "NA"      "0"                    "NA"     "NA"        
## [6,] "NA"            "NA"      "0"                    "NA"     "NA"        
##      MutRatio_Max P_I P_R P  
## [1,] "NA"         "0" "0" "0"
## [2,] "NA"         "0" "0" "0"
## [3,] "NA"         "0" "0" "0"
## [4,] "NA"         "0" "0" "0"
## [5,] "NA"         "0" "0" "0"
## [6,] "NA"         "0" "0" "0"
```


``` r
  print(Export_Summary_IndelSV(Input = Result_HLA2_INDEL, , Mut_Rank_th = 10))
```

```
##              Num_All_Alteration        Num_Evaluated_Alteration 
##                               3                               3 
## Num_Alteration_Generating_NeoAg                 Num_All_Peptide 
##                               2                              46 
##           Num_Evaluated_Peptide    Num_Peptide_Generating_NeoAg 
##                              46                              10
```


``` r
  print(Export_Summary_IndelSV_perFragments(Input = Result_HLA2_INDEL, , Mut_Rank_th = 5))
```

```
##                                                     YEAIYHGIPMVGIPFVLGSTX-0_UGT2B28
## Num_Peptide_Per_Pep                                                           6.000
## Num_Cond_Peptide_Per_Pep                                                      6.000
## Num_Rest_Peptide_Per_Pep                                                      0.000
## Num_Rest_Peptide_Per_Pep / Num_Cond_Peptide_Per_Pep                           0.000
## -logP                                                                         0.483
##                                                     QYYQRDELLAAIKKFHIKVATQX-1_COL12A1
## Num_Peptide_Per_Pep                                                             8.000
## Num_Cond_Peptide_Per_Pep                                                        8.000
## Num_Rest_Peptide_Per_Pep                                                        5.000
## Num_Rest_Peptide_Per_Pep / Num_Cond_Peptide_Per_Pep                             0.625
## -logP                                                                           0.529
##                                                     IMEMARDFLPSLKNRFWKPSILPIFMYKHCSVQFSVRHGDVQTKVHX-2_SLCO1C1
## Num_Peptide_Per_Pep                                                                                    32.000
## Num_Cond_Peptide_Per_Pep                                                                               32.000
## Num_Rest_Peptide_Per_Pep                                                                                4.000
## Num_Rest_Peptide_Per_Pep / Num_Cond_Peptide_Per_Pep                                                     0.125
## -logP                                                                                                   1.081
```

To calculate the binding affinity of neoantigen candaites, which are generated from from SVs, to HLA ClassI.


``` r
  Result_HLA1_SV <- MainSVFUSIONClass1(input_file = "./../data/sample_sv_bnd.txt",
                                       file_name_in_hla_table = "sample",
                                       hla_file = "./../data/sample_hla_table_c1.txt",
                                       refflat_file  = "./../lib/refFlat.grch37.txt",
                                       refmrna_file = "./../lib/refMrna.grch37.fa",
                                       rnaexp_file = "./../data/sample_rna_exp.txt",
                                       netMHCpan_dir = "./../lib/netMHCpan-4.1/netMHCpan",
                                       refdna_file = "./../lib/GRCh37.fa",
                                       mutation_alt_bnd_column = 5,
                                       gene_symbol_column = 7,
                                       mate_id_column = 8)
```

```
## [1] "HLAtype: A*02:01" "HLAtype: A*32:01" "HLAtype: B*15:17" "HLAtype: B*51:01"
## [5] "HLAtype: C*07:01" "HLAtype: C*15:02"
## [1] "Set chr_column as 1 Chr"
## [1] "Set mutation_start_column as 2 Start"
## [1] "Set mutation_end_column as 3 End"
## [1] "Set mutation_ref_column as 4 Ref"
## [1] "Set mutation_alt_column as 5 Alt"
## [1] "Please Manually Indicate depth_normal_column if you want to use."
## [1] "Set depth_normal_column as NA NA"
## [1] "Please Manually Indicate depth_tumor_column if you want to use."
## [1] "Set depth_tumor_column as NA NA"
## [1] "Set nm_id_column as 0 "
## [1] "Please Confirm that Reading Start Line is 1"
## [1] "Chr\tStart\tEnd\tRef\tAlt\tFunc.refGene\tGene.refGene\tmateID"
## [1] "Start Analysis: Mutation SVMERGE137"
## [1] "OK SVMERGE137"
## [1] "Start Analysis: Mutation SVMERGE15"
## [1] "Annotation is Intronic but Position is not Intron"
## [1] "Annotation is Intronic but Position is not Intron"
## [1] "OK SVMERGE15"
## [1] "Start Analysis: Mutation SVMERGE116"
## [1] "OK SVMERGE116"
## [1] "Start Analysis: Mutation SVMERGE3"
## [1] "Annotation is Exonic but Position is Intronic"
## [1] "Annotation is Exonic but Position is Intronic"
## [1] "Annotation is Intronic but Position is not Intron"
## [1] "Annotation is Intronic but Position is not Intron"
## [1] "OK SVMERGE3"
## [1] "Start Analysis: Mutation SVMERGE38"
## [1] "Annotation is Intronic but Position is not Intron"
## [1] "Annotation is Intronic but Position is not Intron"
## [1] "Annotation is Intronic but Position is not Intron"
## [1] "Annotation is Intronic but Position is not Intron"
## [1] "Annotation is Intronic but Position is not Intron"
## [1] "Annotation is Intronic but Position is not Intron"
## [1] "OK SVMERGE38"
## [1] "OK SVMERGE38"
## [1] "OK SVMERGE38"
## [1] "OK SVMERGE38"
## [1] "Start Analysis: Mutation SVMERGE195"
## [1] "OK SVMERGE195"
## [1] "OK SVMERGE195"
## [1] "Start Analysis: Mutation SVMERGE10"
## [1] "Start Analysis: Mutation SVMERGE272"
## [1] "OK SVMERGE272"
## [1] "Start Analysis: Mutation SVMERGE101"
## [1] "Annotation is Intronic but Position is not Intron"
## [1] "Annotation is Intronic but Position is not Intron"
## [1] "Annotation is Intronic but Position is not Intron"
## [1] "Annotation is Intronic but Position is not Intron"
## [1] "Annotation is Intronic but Position is not Intron"
## [1] "Annotation is Intronic but Position is not Intron"
## [1] "Annotation is Intronic but Position is not Intron"
## [1] "Annotation is Intronic but Position is not Intron"
## [1] "OK SVMERGE101"
## [1] "Start Analysis: Mutation SVMERGE40"
## [1] "OK SVMERGE40"
## [1] "Start Analysis: Mutation SVMERGE460"
## [1] "OK SVMERGE460"
## [1] "Start Analysis: Mutation SVMERGE133"
## [1] "Annotation is Exonic but Position is Intronic"
## [1] "Annotation is Exonic but Position is Intronic"
## [1] "OK SVMERGE133"
## [1] "Start Analysis: Mutation SVMERGE31"
## [1] "Annotation is Exonic but Position is Intronic"
## [1] "Annotation is Exonic but Position is Intronic"
## [1] "Start Analysis: Mutation SVMERGE157"
## [1] "Annotation is Exonic but Position is Intronic"
## [1] "Annotation is Exonic but Position is Intronic"
## [1] "Start Analysis: Mutation SVMERGE72"
## [1] "Annotation is Exonic but Position is Intronic"
## [1] "Annotation is Exonic but Position is Intronic"
## [1] "Start Analysis: Mutation SVMERGE177"
## [1] "Annotation is Exonic but Position is Intronic"
## [1] "Annotation is Exonic but Position is Intronic"
## [1] "Start Analysis: Mutation SVMERGE37"
## [1] "Non-Stop!!"
## [1] "Non-Stop!!"
## [1] "OK SVMERGE37"
## [1] "Start Analysis: Mutation SVMERGE206"
## [1] "OK SVMERGE206"
## [1] "Start Analysis: Mutation SVMERGE58"
## [1] "OK SVMERGE58"
## [1] "Start Analysis: Mutation SVMERGE210"
## [1] "OK SVMERGE210"
## [1] "Start Analysis: Mutation SVMERGE205"
## [1] "Annotation is Intronic but Position is not Intron"
## [1] "Annotation is Intronic but Position is not Intron"
## [1] "OK SVMERGE205"
## [1] "OK SVMERGE205"
## [1] "OK SVMERGE205"
## [1] "OK SVMERGE205"
## [1] "OK SVMERGE205"
## [1] "result.ID.SV1/sample_sv_bnd.txt.ID_SVFusion.peptide.txt.list.vcf"
## [1] "Executing netMHCpan to result.ID.SV1"
## [1] "./../lib/netMHCpan-4.1/netMHCpan -BA  -l 8,9,10,11,12,13 -f result.ID.SV1/sample_sv_bnd.txt.ID_SVFusion.peptide.fasta -a HLA-A02:01 > result.ID.SV1/sample_sv_bnd.txt.ID_SVFusion.HLACLASS1.1.peptide.txt"
## [1] "./../lib/netMHCpan-4.1/netMHCpan -BA  -l 8,9,10,11,12,13 -f result.ID.SV1/sample_sv_bnd.txt.ID_SVFusion.peptide.fasta -a HLA-A32:01 > result.ID.SV1/sample_sv_bnd.txt.ID_SVFusion.HLACLASS1.2.peptide.txt"
## [1] "./../lib/netMHCpan-4.1/netMHCpan -BA  -l 8,9,10,11,12,13 -f result.ID.SV1/sample_sv_bnd.txt.ID_SVFusion.peptide.fasta -a HLA-B15:17 > result.ID.SV1/sample_sv_bnd.txt.ID_SVFusion.HLACLASS1.3.peptide.txt"
## [1] "./../lib/netMHCpan-4.1/netMHCpan -BA  -l 8,9,10,11,12,13 -f result.ID.SV1/sample_sv_bnd.txt.ID_SVFusion.peptide.fasta -a HLA-B51:01 > result.ID.SV1/sample_sv_bnd.txt.ID_SVFusion.HLACLASS1.4.peptide.txt"
## [1] "./../lib/netMHCpan-4.1/netMHCpan -BA  -l 8,9,10,11,12,13 -f result.ID.SV1/sample_sv_bnd.txt.ID_SVFusion.peptide.fasta -a HLA-C07:01 > result.ID.SV1/sample_sv_bnd.txt.ID_SVFusion.HLACLASS1.5.peptide.txt"
## [1] "./../lib/netMHCpan-4.1/netMHCpan -BA  -l 8,9,10,11,12,13 -f result.ID.SV1/sample_sv_bnd.txt.ID_SVFusion.peptide.fasta -a HLA-C15:02 > result.ID.SV1/sample_sv_bnd.txt.ID_SVFusion.HLACLASS1.6.peptide.txt"
## [1] "Merging Results..."
## [1] "/Users/takanorihasegawa/Documents/GitHub/Neoantimon/vignettes/result.ID.SV1/sample_sv_bnd.txt.ID_SVFusion.HLACLASS1.1.peptide.txt"
## [1] "4.16666666666667 perc. fin"
## [1] "8.33333333333333 perc. fin"
## [1] "12.5 perc. fin"
## [1] "16.6666666666667 perc. fin"
## [1] "20.8333333333333 perc. fin"
## [1] "25 perc. fin"
## [1] "29.1666666666667 perc. fin"
## [1] "33.3333333333333 perc. fin"
## [1] "37.5 perc. fin"
## [1] "41.6666666666667 perc. fin"
## [1] "45.8333333333333 perc. fin"
## [1] "50 perc. fin"
## [1] "54.1666666666667 perc. fin"
## [1] "58.3333333333333 perc. fin"
## [1] "62.5 perc. fin"
## [1] "66.6666666666667 perc. fin"
## [1] "70.8333333333333 perc. fin"
## [1] "75 perc. fin"
## [1] "79.1666666666667 perc. fin"
## [1] "83.3333333333333 perc. fin"
## [1] "87.5 perc. fin"
## [1] "91.6666666666667 perc. fin"
## [1] "95.8333333333333 perc. fin"
## [1] "100 perc. fin"
## [1] "/Users/takanorihasegawa/Documents/GitHub/Neoantimon/vignettes/result.ID.SV1/sample_sv_bnd.txt.ID_SVFusion.HLACLASS1.2.peptide.txt"
## [1] "4.16666666666667 perc. fin"
## [1] "8.33333333333333 perc. fin"
## [1] "12.5 perc. fin"
## [1] "16.6666666666667 perc. fin"
## [1] "20.8333333333333 perc. fin"
## [1] "25 perc. fin"
## [1] "29.1666666666667 perc. fin"
## [1] "33.3333333333333 perc. fin"
## [1] "37.5 perc. fin"
## [1] "41.6666666666667 perc. fin"
## [1] "45.8333333333333 perc. fin"
## [1] "50 perc. fin"
## [1] "54.1666666666667 perc. fin"
## [1] "58.3333333333333 perc. fin"
## [1] "62.5 perc. fin"
## [1] "66.6666666666667 perc. fin"
## [1] "70.8333333333333 perc. fin"
## [1] "75 perc. fin"
## [1] "79.1666666666667 perc. fin"
## [1] "83.3333333333333 perc. fin"
## [1] "87.5 perc. fin"
## [1] "91.6666666666667 perc. fin"
## [1] "95.8333333333333 perc. fin"
## [1] "100 perc. fin"
## [1] "/Users/takanorihasegawa/Documents/GitHub/Neoantimon/vignettes/result.ID.SV1/sample_sv_bnd.txt.ID_SVFusion.HLACLASS1.3.peptide.txt"
## [1] "4.16666666666667 perc. fin"
## [1] "8.33333333333333 perc. fin"
## [1] "12.5 perc. fin"
## [1] "16.6666666666667 perc. fin"
## [1] "20.8333333333333 perc. fin"
## [1] "25 perc. fin"
## [1] "29.1666666666667 perc. fin"
## [1] "33.3333333333333 perc. fin"
## [1] "37.5 perc. fin"
## [1] "41.6666666666667 perc. fin"
## [1] "45.8333333333333 perc. fin"
## [1] "50 perc. fin"
## [1] "54.1666666666667 perc. fin"
## [1] "58.3333333333333 perc. fin"
## [1] "62.5 perc. fin"
## [1] "66.6666666666667 perc. fin"
## [1] "70.8333333333333 perc. fin"
## [1] "75 perc. fin"
## [1] "79.1666666666667 perc. fin"
## [1] "83.3333333333333 perc. fin"
## [1] "87.5 perc. fin"
## [1] "91.6666666666667 perc. fin"
## [1] "95.8333333333333 perc. fin"
## [1] "100 perc. fin"
## [1] "/Users/takanorihasegawa/Documents/GitHub/Neoantimon/vignettes/result.ID.SV1/sample_sv_bnd.txt.ID_SVFusion.HLACLASS1.4.peptide.txt"
## [1] "4.16666666666667 perc. fin"
## [1] "8.33333333333333 perc. fin"
## [1] "12.5 perc. fin"
## [1] "16.6666666666667 perc. fin"
## [1] "20.8333333333333 perc. fin"
## [1] "25 perc. fin"
## [1] "29.1666666666667 perc. fin"
## [1] "33.3333333333333 perc. fin"
## [1] "37.5 perc. fin"
## [1] "41.6666666666667 perc. fin"
## [1] "45.8333333333333 perc. fin"
## [1] "50 perc. fin"
## [1] "54.1666666666667 perc. fin"
## [1] "58.3333333333333 perc. fin"
## [1] "62.5 perc. fin"
## [1] "66.6666666666667 perc. fin"
## [1] "70.8333333333333 perc. fin"
## [1] "75 perc. fin"
## [1] "79.1666666666667 perc. fin"
## [1] "83.3333333333333 perc. fin"
## [1] "87.5 perc. fin"
## [1] "91.6666666666667 perc. fin"
## [1] "95.8333333333333 perc. fin"
## [1] "100 perc. fin"
## [1] "/Users/takanorihasegawa/Documents/GitHub/Neoantimon/vignettes/result.ID.SV1/sample_sv_bnd.txt.ID_SVFusion.HLACLASS1.5.peptide.txt"
## [1] "4.16666666666667 perc. fin"
## [1] "8.33333333333333 perc. fin"
## [1] "12.5 perc. fin"
## [1] "16.6666666666667 perc. fin"
## [1] "20.8333333333333 perc. fin"
## [1] "25 perc. fin"
## [1] "29.1666666666667 perc. fin"
## [1] "33.3333333333333 perc. fin"
## [1] "37.5 perc. fin"
## [1] "41.6666666666667 perc. fin"
## [1] "45.8333333333333 perc. fin"
## [1] "50 perc. fin"
## [1] "54.1666666666667 perc. fin"
## [1] "58.3333333333333 perc. fin"
## [1] "62.5 perc. fin"
## [1] "66.6666666666667 perc. fin"
## [1] "70.8333333333333 perc. fin"
## [1] "75 perc. fin"
## [1] "79.1666666666667 perc. fin"
## [1] "83.3333333333333 perc. fin"
## [1] "87.5 perc. fin"
## [1] "91.6666666666667 perc. fin"
## [1] "95.8333333333333 perc. fin"
## [1] "100 perc. fin"
## [1] "/Users/takanorihasegawa/Documents/GitHub/Neoantimon/vignettes/result.ID.SV1/sample_sv_bnd.txt.ID_SVFusion.HLACLASS1.6.peptide.txt"
## [1] "4.16666666666667 perc. fin"
## [1] "8.33333333333333 perc. fin"
## [1] "12.5 perc. fin"
## [1] "16.6666666666667 perc. fin"
## [1] "20.8333333333333 perc. fin"
## [1] "25 perc. fin"
## [1] "29.1666666666667 perc. fin"
## [1] "33.3333333333333 perc. fin"
## [1] "37.5 perc. fin"
## [1] "41.6666666666667 perc. fin"
## [1] "45.8333333333333 perc. fin"
## [1] "50 perc. fin"
## [1] "54.1666666666667 perc. fin"
## [1] "58.3333333333333 perc. fin"
## [1] "62.5 perc. fin"
## [1] "66.6666666666667 perc. fin"
## [1] "70.8333333333333 perc. fin"
## [1] "75 perc. fin"
## [1] "79.1666666666667 perc. fin"
## [1] "83.3333333333333 perc. fin"
## [1] "87.5 perc. fin"
## [1] "91.6666666666667 perc. fin"
## [1] "95.8333333333333 perc. fin"
## [1] "100 perc. fin"
## [1] "Successfully Finished."
```

``` r
  Result_HLA1_SV <- CalculatePriorityScores(result = Result_HLA1_SV, useRNAvaf = FALSE)
  print(head(Result_HLA1_SV))
```

```
##      HLA           Pos Gene     Evaluated_Mutant_Peptide_Core
## [1,] "HLA-A*02:01" "1" "0_AAR2" "DYNSWEVGPKFRE"              
## [2,] "HLA-A*02:01" "2" "0_AAR2" "YNSWEVGPKFREQ"              
## [3,] "HLA-A*02:01" "3" "0_AAR2" "NSWEVGPKFREQL"              
## [4,] "HLA-A*02:01" "4" "0_AAR2" "SWEVGPKFREQLK"              
## [5,] "HLA-A*02:01" "5" "0_AAR2" "WEVGPKFREQLKL"              
## [6,] "HLA-A*02:01" "6" "0_AAR2" "EVGPKFREQLKL"               
##      Evaluated_Mutant_Peptide Mut_EL      Mut_Rank Chr  NM_ID                
## [1,] "DYNSWEVGPKFRE"          "0.0000010" "95.000" "20" "NM_015511_NM_015906"
## [2,] "YNSWEVGPKFREQ"          "0.0000080" "76.667" "20" "NM_015511_NM_015906"
## [3,] "NSWEVGPKFREQL"          "0.0004280" "24.750" "20" "NM_015511_NM_015906"
## [4,] "SWEVGPKFREQLK"          "0.0000030" "87.500" "20" "NM_015511_NM_015906"
## [5,] "WEVGPKFREQLKL"          "0.0000430" "53.333" "20" "NM_015511_NM_015906"
## [6,] "EVGPKFREQLKLF"          "0.0000060" "80.000" "20" "NM_015511_NM_015906"
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
##      P_I                 P_R                     P                      
## [1,] "0.994835363484567" "1.12475327099476e-202" "1.12475613722433e-202"
## [2,] "0.994835363426825" "7.25539311432306e-163" "7.25541160337557e-163"
## [3,] "0.994835359962294" "3.95361188289843e-50"  "3.95362195796073e-50" 
## [4,] "0.994835363468069" "2.17320320184772e-186" "2.17320873986164e-186"
## [5,] "0.994835363138115" "3.38686707516126e-112" "3.38687570597729e-112"
## [6,] "0.994835363443322" "4.19897614731469e-170" "4.19898684764329e-170"
```


``` r
  print(Export_Summary_IndelSV(Result_HLA1_SV, , Mut_Rank_th = 5))
```

```
##              Num_All_Alteration        Num_Evaluated_Alteration 
##                              16                              16 
## Num_Alteration_Generating_NeoAg                 Num_All_Peptide 
##                              16                            3183 
##           Num_Evaluated_Peptide    Num_Peptide_Generating_NeoAg 
##                            3183                             351
```


``` r
  print(Export_Summary_IndelSV_perFragments(Result_HLA1_SV, , Mut_Rank_th = 5))
```

```
##                                                     DYNSWEVGPKFREQLKLFCETCDRLTCRDCQLLE-0_AAR2
## Num_Peptide_Per_Pep                                                                    57.000
## Num_Cond_Peptide_Per_Pep                                                               57.000
## Num_Rest_Peptide_Per_Pep                                                                5.000
## Num_Rest_Peptide_Per_Pep / Num_Cond_Peptide_Per_Pep                                     0.088
## -logP                                                                                   0.782
##                                                     GDPHELRNIFLQLSAVQEAQLKRLEVTRPRVLGSREQGQVPRMARQPPPPWVHAAFLLCLLSLGG-1_SLC25A12
## Num_Peptide_Per_Pep                                                                                                      237.000
## Num_Cond_Peptide_Per_Pep                                                                                                 237.000
## Num_Rest_Peptide_Per_Pep                                                                                                  31.000
## Num_Rest_Peptide_Per_Pep / Num_Cond_Peptide_Per_Pep                                                                        0.131
## -logP                                                                                                                      1.496
##                                                     TMAEKRQLFIEMLYRLX-2_DTNB
## Num_Peptide_Per_Pep                                                   24.000
## Num_Cond_Peptide_Per_Pep                                              24.000
## Num_Rest_Peptide_Per_Pep                                               7.000
## Num_Rest_Peptide_Per_Pep / Num_Cond_Peptide_Per_Pep                    0.292
## -logP                                                                  0.391
##                                                     TDWLSDCCFHPSERCRCTLYGHTDSVNSIEFFP-3_SPAG16
## Num_Peptide_Per_Pep                                                                     51.000
## Num_Cond_Peptide_Per_Pep                                                                51.000
## Num_Rest_Peptide_Per_Pep                                                                 2.000
## Num_Rest_Peptide_Per_Pep / Num_Cond_Peptide_Per_Pep                                      0.039
## -logP                                                                                    0.759
##                                                     MKSCKRNLSTTTGSSERX-4_TMCC1
## Num_Peptide_Per_Pep                                                     45.000
## Num_Cond_Peptide_Per_Pep                                                45.000
## Num_Rest_Peptide_Per_Pep                                                 0.000
## Num_Rest_Peptide_Per_Pep / Num_Cond_Peptide_Per_Pep                      0.000
## -logP                                                                    0.414
##                                                     ACLPGEEGTAERSCKRNLSTTTGSSERX-5_TMCC1
## Num_Peptide_Per_Pep                                                               90.000
## Num_Cond_Peptide_Per_Pep                                                          90.000
## Num_Rest_Peptide_Per_Pep                                                           0.000
## Num_Rest_Peptide_Per_Pep / Num_Cond_Peptide_Per_Pep                                0.000
## -logP                                                                              0.644
##                                                     VQRFSLRRQLSKSCKRNLSTTTGSSERX-6_TMCC1
## Num_Peptide_Per_Pep                                                               90.000
## Num_Cond_Peptide_Per_Pep                                                          90.000
## Num_Rest_Peptide_Per_Pep                                                           3.000
## Num_Rest_Peptide_Per_Pep / Num_Cond_Peptide_Per_Pep                                0.033
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
## Num_Rest_Peptide_Per_Pep                                                                12.000
## Num_Rest_Peptide_Per_Pep / Num_Cond_Peptide_Per_Pep                                      0.100
## -logP                                                                                    0.759
##                                                     GQLRGLIPSSPQKTLECKQQKIILLRREMKETX-9_PDLIM3
## Num_Peptide_Per_Pep                                                                    120.000
## Num_Cond_Peptide_Per_Pep                                                               120.000
## Num_Rest_Peptide_Per_Pep                                                                19.000
## Num_Rest_Peptide_Per_Pep / Num_Cond_Peptide_Per_Pep                                      0.158
## -logP                                                                                    0.759
##                                                     GDQYIVDMANTKRTLLWTRPMTVILGKPFCVMPKQQKTAHIGFLQHIPRLSPKPCLPKLNLMMRKQRMSQNGKNVKFEESHLRAVCMSGRGMGQVWVFFLCSX-10_WDR70
## Num_Peptide_Per_Pep                                                                                                                                           540.00
## Num_Cond_Peptide_Per_Pep                                                                                                                                      540.00
## Num_Rest_Peptide_Per_Pep                                                                                                                                       81.00
## Num_Rest_Peptide_Per_Pep / Num_Cond_Peptide_Per_Pep                                                                                                             0.15
## -logP                                                                                                                                                           2.37
##                                                     LCPASRALEEKKGHELLGGPSLGAPRPGSQERTGENTAACQDHRFWAGQTAGCGRERIPCRRRQSAYQVDGIGINFTQNLYPPEX-11_EGFR
## Num_Peptide_Per_Pep                                                                                                                       432.000
## Num_Cond_Peptide_Per_Pep                                                                                                                  432.000
## Num_Rest_Peptide_Per_Pep                                                                                                                   24.000
## Num_Rest_Peptide_Per_Pep / Num_Cond_Peptide_Per_Pep                                                                                         0.056
## -logP                                                                                                                                       1.956
##                                                     LTPMPEGLSQQQDMLKNTSKGHPDRLPLQMALT-12_ARHGEF10
## Num_Peptide_Per_Pep                                                                        51.000
## Num_Cond_Peptide_Per_Pep                                                                   51.000
## Num_Rest_Peptide_Per_Pep                                                                    3.000
## Num_Rest_Peptide_Per_Pep / Num_Cond_Peptide_Per_Pep                                         0.059
## -logP                                                                                       0.759
##                                                     PGVQGKKVRKPSPTTSASSSSWTSAGTSWPAGRRTGTATSGTATTTSVWPGCGTRMWSTQWSSVPRSRSCCSRPATTPPSKPGAPHAPCASSRHLAHGLAPSSPGLPARGAEVCWVHWSHRDPLRTSPGSVAFSRAGEVEMLIAVTPX-13_NOTCH1
## Num_Peptide_Per_Pep                                                                                                                                                                                        810.000
## Num_Cond_Peptide_Per_Pep                                                                                                                                                                                   810.000
## Num_Rest_Peptide_Per_Pep                                                                                                                                                                                    97.000
## Num_Rest_Peptide_Per_Pep / Num_Cond_Peptide_Per_Pep                                                                                                                                                          0.120
## -logP                                                                                                                                                                                                        3.405
##                                                     SRGRGNNNRKGREVTPL-14_UBAP2
## Num_Peptide_Per_Pep                                                     30.000
## Num_Cond_Peptide_Per_Pep                                                30.000
## Num_Rest_Peptide_Per_Pep                                                 1.000
## Num_Rest_Peptide_Per_Pep / Num_Cond_Peptide_Per_Pep                      0.033
## -logP                                                                    0.391
##                                                     IPSSQLAAKLLHMLTMRMLSKSATGRWX-15_RASGRP2
## Num_Peptide_Per_Pep                                                                  90.000
## Num_Cond_Peptide_Per_Pep                                                             90.000
## Num_Rest_Peptide_Per_Pep                                                             20.000
## Num_Rest_Peptide_Per_Pep / Num_Cond_Peptide_Per_Pep                                   0.222
## -logP                                                                                 0.644
##                                                     RKSPPEKKLRRYPPGQGATIGVDFMIKTVEINGE-16_GAB2
## Num_Peptide_Per_Pep                                                                     57.000
## Num_Cond_Peptide_Per_Pep                                                                57.000
## Num_Rest_Peptide_Per_Pep                                                                11.000
## Num_Rest_Peptide_Per_Pep / Num_Cond_Peptide_Per_Pep                                      0.193
## -logP                                                                                    0.782
##                                                     MSVIVLPQPDEVLNLVQSYVTLRVPLYVSYVFH-17_MTRF1
## Num_Peptide_Per_Pep                                                                     51.000
## Num_Cond_Peptide_Per_Pep                                                                51.000
## Num_Rest_Peptide_Per_Pep                                                                12.000
## Num_Rest_Peptide_Per_Pep / Num_Cond_Peptide_Per_Pep                                      0.235
## -logP                                                                                    0.759
##                                                     LWKCLRKPVLNDRNLQLHTDKGSFLKEKNKKLKKK-18_EVI2B
## Num_Peptide_Per_Pep                                                                       63.000
## Num_Cond_Peptide_Per_Pep                                                                  63.000
## Num_Rest_Peptide_Per_Pep                                                                   8.000
## Num_Rest_Peptide_Per_Pep / Num_Cond_Peptide_Per_Pep                                        0.127
## -logP                                                                                      0.805
##                                                     RTFLQSTDPGDTGGVFEGSTLCPSDPAFSRTDDPX-19_IKZF3
## Num_Peptide_Per_Pep                                                                      132.000
## Num_Cond_Peptide_Per_Pep                                                                 132.000
## Num_Rest_Peptide_Per_Pep                                                                   9.000
## Num_Rest_Peptide_Per_Pep / Num_Cond_Peptide_Per_Pep                                        0.068
## -logP                                                                                      0.805
##                                                     RDALTGHLRTHSGGVFEGSTLCPSDPAFSRTDDPX-20_IKZF3
## Num_Peptide_Per_Pep                                                                      132.000
## Num_Cond_Peptide_Per_Pep                                                                 132.000
## Num_Rest_Peptide_Per_Pep                                                                  10.000
## Num_Rest_Peptide_Per_Pep / Num_Cond_Peptide_Per_Pep                                        0.076
## -logP                                                                                      0.805
##                                                     GEGPANEDEDIGGGVFEGSTLCPSDPAFSRTDDPX-21_IKZF3
## Num_Peptide_Per_Pep                                                                      132.000
## Num_Cond_Peptide_Per_Pep                                                                 132.000
## Num_Rest_Peptide_Per_Pep                                                                   5.000
## Num_Rest_Peptide_Per_Pep / Num_Cond_Peptide_Per_Pep                                        0.038
## -logP                                                                                      0.805
##                                                     NVLMVHKRSHTGGVFEGSTLCPSDPAFSRTDDPX-22_IKZF3
## Num_Peptide_Per_Pep                                                                     126.000
## Num_Cond_Peptide_Per_Pep                                                                126.000
## Num_Rest_Peptide_Per_Pep                                                                 10.000
## Num_Rest_Peptide_Per_Pep / Num_Cond_Peptide_Per_Pep                                       0.079
## -logP                                                                                     0.782
##                                                     FNVLMVHKRSHTGGVFEGSTLCPSDPAFSRTDDPX-23_IKZF3
## Num_Peptide_Per_Pep                                                                      132.000
## Num_Cond_Peptide_Per_Pep                                                                 132.000
## Num_Rest_Peptide_Per_Pep                                                                  10.000
## Num_Rest_Peptide_Per_Pep / Num_Cond_Peptide_Per_Pep                                        0.076
## -logP                                                                                      0.805
```

To calculate the binding affinity of neoantigen candaites, which are generated from from SVs, to HLA ClassII.


``` r
  Result_HLA2_SV <- MainSVFUSIONClass2(input_file = "./../data/sample_sv_bnd.txt",
                                       file_name_in_hla_table = "sample",
                                       hla_file = "./../data/sample_hla_table_c2.txt",
                                       refflat_file  = "./../lib/refFlat.grch37.txt",
                                       refmrna_file = "./../lib/refMrna.grch37.fa",
                                       rnaexp_file = "./../data/sample_rna_exp.txt",
                                       netMHCIIpan_dir = "./../lib/netMHCIIpan-4.3/netMHCIIpan",
                                       refdna_file = "./../lib/GRCh37.fa",
                                       mutation_alt_bnd_column = 5,
                                       gene_symbol_column = 7,
                                       mate_id_column = 8)
```

```
##  [1] "HLAtype: DPA1*01:03" "HLAtype: DPA1*02:01" "HLAtype: DPB1*02:01"
##  [4] "HLAtype: DPB1*09:01" "HLAtype: DQA1*01:02" "HLAtype: DQA1*05:05"
##  [7] "HLAtype: DQB1*03:01" "HLAtype: DQB1*06:04" "HLAtype: DRB1*11:04"
## [10] "HLAtype: DRB1*13:02"
## [1] "Set chr_column as 1 Chr"
## [1] "Set mutation_start_column as 2 Start"
## [1] "Set mutation_end_column as 3 End"
## [1] "Set mutation_ref_column as 4 Ref"
## [1] "Set mutation_alt_column as 5 Alt"
## [1] "Please Manually Indicate depth_normal_column if you want to use."
## [1] "Set depth_normal_column as NA NA"
## [1] "Please Manually Indicate depth_tumor_column if you want to use."
## [1] "Set depth_tumor_column as NA NA"
## [1] "Set nm_id_column as 0 "
## [1] "Please Confirm that Reading Start Line is 1"
## [1] "Chr\tStart\tEnd\tRef\tAlt\tFunc.refGene\tGene.refGene\tmateID"
## [1] "Start Analysis: Mutation SVMERGE137"
## [1] "OK SVMERGE137"
## [1] "Start Analysis: Mutation SVMERGE15"
## [1] "Annotation is Intronic but Position is not Intron"
## [1] "Annotation is Intronic but Position is not Intron"
## [1] "OK SVMERGE15"
## [1] "Start Analysis: Mutation SVMERGE116"
## [1] "OK SVMERGE116"
## [1] "Start Analysis: Mutation SVMERGE3"
## [1] "Annotation is Exonic but Position is Intronic"
## [1] "Annotation is Exonic but Position is Intronic"
## [1] "Annotation is Intronic but Position is not Intron"
## [1] "Annotation is Intronic but Position is not Intron"
## [1] "OK SVMERGE3"
## [1] "Start Analysis: Mutation SVMERGE38"
## [1] "Annotation is Intronic but Position is not Intron"
## [1] "Annotation is Intronic but Position is not Intron"
## [1] "Annotation is Intronic but Position is not Intron"
## [1] "Annotation is Intronic but Position is not Intron"
## [1] "Annotation is Intronic but Position is not Intron"
## [1] "Annotation is Intronic but Position is not Intron"
## [1] "OK SVMERGE38"
## [1] "OK SVMERGE38"
## [1] "OK SVMERGE38"
## [1] "OK SVMERGE38"
## [1] "Start Analysis: Mutation SVMERGE195"
## [1] "OK SVMERGE195"
## [1] "OK SVMERGE195"
## [1] "Start Analysis: Mutation SVMERGE10"
## [1] "Start Analysis: Mutation SVMERGE272"
## [1] "OK SVMERGE272"
## [1] "Start Analysis: Mutation SVMERGE101"
## [1] "Annotation is Intronic but Position is not Intron"
## [1] "Annotation is Intronic but Position is not Intron"
## [1] "Annotation is Intronic but Position is not Intron"
## [1] "Annotation is Intronic but Position is not Intron"
## [1] "Annotation is Intronic but Position is not Intron"
## [1] "Annotation is Intronic but Position is not Intron"
## [1] "Annotation is Intronic but Position is not Intron"
## [1] "Annotation is Intronic but Position is not Intron"
## [1] "OK SVMERGE101"
## [1] "Start Analysis: Mutation SVMERGE40"
## [1] "OK SVMERGE40"
## [1] "Start Analysis: Mutation SVMERGE460"
## [1] "OK SVMERGE460"
## [1] "Start Analysis: Mutation SVMERGE133"
## [1] "Annotation is Exonic but Position is Intronic"
## [1] "Annotation is Exonic but Position is Intronic"
## [1] "OK SVMERGE133"
## [1] "Start Analysis: Mutation SVMERGE31"
## [1] "Annotation is Exonic but Position is Intronic"
## [1] "Annotation is Exonic but Position is Intronic"
## [1] "Start Analysis: Mutation SVMERGE157"
## [1] "Annotation is Exonic but Position is Intronic"
## [1] "Annotation is Exonic but Position is Intronic"
## [1] "Start Analysis: Mutation SVMERGE72"
## [1] "Annotation is Exonic but Position is Intronic"
## [1] "Annotation is Exonic but Position is Intronic"
## [1] "Start Analysis: Mutation SVMERGE177"
## [1] "Annotation is Exonic but Position is Intronic"
## [1] "Annotation is Exonic but Position is Intronic"
## [1] "Start Analysis: Mutation SVMERGE37"
## [1] "Non-Stop!!"
## [1] "Non-Stop!!"
## [1] "OK SVMERGE37"
## [1] "Start Analysis: Mutation SVMERGE206"
## [1] "OK SVMERGE206"
## [1] "Start Analysis: Mutation SVMERGE58"
## [1] "OK SVMERGE58"
## [1] "Start Analysis: Mutation SVMERGE210"
## [1] "OK SVMERGE210"
## [1] "Start Analysis: Mutation SVMERGE205"
## [1] "Annotation is Intronic but Position is not Intron"
## [1] "Annotation is Intronic but Position is not Intron"
## [1] "OK SVMERGE205"
## [1] "OK SVMERGE205"
## [1] "OK SVMERGE205"
## [1] "OK SVMERGE205"
## [1] "OK SVMERGE205"
## [1] "result.ID.SV2/sample_sv_bnd.txt.ID_SVFusion.peptide.txt.list.vcf"
## [1] "Executing netMHCpan to result.ID.SV2"
## [1] "./../lib/netMHCIIpan-4.3/netMHCIIpan -length 15 -f result.ID.SV2/sample_sv_bnd.txt.ID_SVFusion.peptide.fasta -choose -cha DPA10103 -choose -chb DPB10201 > result.ID.SV2/sample_sv_bnd.txt.ID_SVFusion.HLACLASS2.1.peptide.txt"
## [1] "./../lib/netMHCIIpan-4.3/netMHCIIpan -length 15 -f result.ID.SV2/sample_sv_bnd.txt.ID_SVFusion.peptide.fasta -choose -cha DPA10103 -choose -chb DPB10901 > result.ID.SV2/sample_sv_bnd.txt.ID_SVFusion.HLACLASS2.2.peptide.txt"
## [1] "./../lib/netMHCIIpan-4.3/netMHCIIpan -length 15 -f result.ID.SV2/sample_sv_bnd.txt.ID_SVFusion.peptide.fasta -choose -cha DPA10201 -choose -chb DPB10201 > result.ID.SV2/sample_sv_bnd.txt.ID_SVFusion.HLACLASS2.3.peptide.txt"
## [1] "./../lib/netMHCIIpan-4.3/netMHCIIpan -length 15 -f result.ID.SV2/sample_sv_bnd.txt.ID_SVFusion.peptide.fasta -choose -cha DPA10201 -choose -chb DPB10901 > result.ID.SV2/sample_sv_bnd.txt.ID_SVFusion.HLACLASS2.4.peptide.txt"
## [1] "./../lib/netMHCIIpan-4.3/netMHCIIpan -length 15 -f result.ID.SV2/sample_sv_bnd.txt.ID_SVFusion.peptide.fasta -choose -cha DQA10102 -choose -chb DQB10301 > result.ID.SV2/sample_sv_bnd.txt.ID_SVFusion.HLACLASS2.5.peptide.txt"
## [1] "./../lib/netMHCIIpan-4.3/netMHCIIpan -length 15 -f result.ID.SV2/sample_sv_bnd.txt.ID_SVFusion.peptide.fasta -choose -cha DQA10102 -choose -chb DQB10604 > result.ID.SV2/sample_sv_bnd.txt.ID_SVFusion.HLACLASS2.6.peptide.txt"
## [1] "./../lib/netMHCIIpan-4.3/netMHCIIpan -length 15 -f result.ID.SV2/sample_sv_bnd.txt.ID_SVFusion.peptide.fasta -choose -cha DQA10505 -choose -chb DQB10301 > result.ID.SV2/sample_sv_bnd.txt.ID_SVFusion.HLACLASS2.7.peptide.txt"
## [1] "./../lib/netMHCIIpan-4.3/netMHCIIpan -length 15 -f result.ID.SV2/sample_sv_bnd.txt.ID_SVFusion.peptide.fasta -choose -cha DQA10505 -choose -chb DQB10604 > result.ID.SV2/sample_sv_bnd.txt.ID_SVFusion.HLACLASS2.8.peptide.txt"
## [1] "./../lib/netMHCIIpan-4.3/netMHCIIpan -length 15 -f result.ID.SV2/sample_sv_bnd.txt.ID_SVFusion.peptide.fasta -a DRB1_1104 > result.ID.SV2/sample_sv_bnd.txt.ID_SVFusion.HLACLASS2.9.peptide.txt"
## [1] "./../lib/netMHCIIpan-4.3/netMHCIIpan -length 15 -f result.ID.SV2/sample_sv_bnd.txt.ID_SVFusion.peptide.fasta -a DRB1_1302 > result.ID.SV2/sample_sv_bnd.txt.ID_SVFusion.HLACLASS2.10.peptide.txt"
## [1] "Merging Results..."
## [1] "/Users/takanorihasegawa/Documents/GitHub/Neoantimon/vignettes/result.ID.SV2/sample_sv_bnd.txt.ID_SVFusion.HLACLASS2.1.peptide.txt"
## [1] "4.16666666666667 perc. fin"
## [1] "8.33333333333333 perc. fin"
## [1] "12.5 perc. fin"
## [1] "16.6666666666667 perc. fin"
## [1] "20.8333333333333 perc. fin"
## [1] "25 perc. fin"
## [1] "29.1666666666667 perc. fin"
## [1] "33.3333333333333 perc. fin"
## [1] "37.5 perc. fin"
## [1] "41.6666666666667 perc. fin"
## [1] "45.8333333333333 perc. fin"
## [1] "50 perc. fin"
## [1] "54.1666666666667 perc. fin"
## [1] "58.3333333333333 perc. fin"
## [1] "62.5 perc. fin"
## [1] "66.6666666666667 perc. fin"
## [1] "70.8333333333333 perc. fin"
## [1] "75 perc. fin"
## [1] "79.1666666666667 perc. fin"
## [1] "83.3333333333333 perc. fin"
## [1] "87.5 perc. fin"
## [1] "91.6666666666667 perc. fin"
## [1] "95.8333333333333 perc. fin"
## [1] "100 perc. fin"
## [1] "/Users/takanorihasegawa/Documents/GitHub/Neoantimon/vignettes/result.ID.SV2/sample_sv_bnd.txt.ID_SVFusion.HLACLASS2.10.peptide.txt"
## [1] "4.16666666666667 perc. fin"
## [1] "8.33333333333333 perc. fin"
## [1] "12.5 perc. fin"
## [1] "16.6666666666667 perc. fin"
## [1] "20.8333333333333 perc. fin"
## [1] "25 perc. fin"
## [1] "29.1666666666667 perc. fin"
## [1] "33.3333333333333 perc. fin"
## [1] "37.5 perc. fin"
## [1] "41.6666666666667 perc. fin"
## [1] "45.8333333333333 perc. fin"
## [1] "50 perc. fin"
## [1] "54.1666666666667 perc. fin"
## [1] "58.3333333333333 perc. fin"
## [1] "62.5 perc. fin"
## [1] "66.6666666666667 perc. fin"
## [1] "70.8333333333333 perc. fin"
## [1] "75 perc. fin"
## [1] "79.1666666666667 perc. fin"
## [1] "83.3333333333333 perc. fin"
## [1] "87.5 perc. fin"
## [1] "91.6666666666667 perc. fin"
## [1] "95.8333333333333 perc. fin"
## [1] "100 perc. fin"
## [1] "/Users/takanorihasegawa/Documents/GitHub/Neoantimon/vignettes/result.ID.SV2/sample_sv_bnd.txt.ID_SVFusion.HLACLASS2.2.peptide.txt"
## [1] "4.16666666666667 perc. fin"
## [1] "8.33333333333333 perc. fin"
## [1] "12.5 perc. fin"
## [1] "16.6666666666667 perc. fin"
## [1] "20.8333333333333 perc. fin"
## [1] "25 perc. fin"
## [1] "29.1666666666667 perc. fin"
## [1] "33.3333333333333 perc. fin"
## [1] "37.5 perc. fin"
## [1] "41.6666666666667 perc. fin"
## [1] "45.8333333333333 perc. fin"
## [1] "50 perc. fin"
## [1] "54.1666666666667 perc. fin"
## [1] "58.3333333333333 perc. fin"
## [1] "62.5 perc. fin"
## [1] "66.6666666666667 perc. fin"
## [1] "70.8333333333333 perc. fin"
## [1] "75 perc. fin"
## [1] "79.1666666666667 perc. fin"
## [1] "83.3333333333333 perc. fin"
## [1] "87.5 perc. fin"
## [1] "91.6666666666667 perc. fin"
## [1] "95.8333333333333 perc. fin"
## [1] "100 perc. fin"
## [1] "/Users/takanorihasegawa/Documents/GitHub/Neoantimon/vignettes/result.ID.SV2/sample_sv_bnd.txt.ID_SVFusion.HLACLASS2.3.peptide.txt"
## [1] "4.16666666666667 perc. fin"
## [1] "8.33333333333333 perc. fin"
## [1] "12.5 perc. fin"
## [1] "16.6666666666667 perc. fin"
## [1] "20.8333333333333 perc. fin"
## [1] "25 perc. fin"
## [1] "29.1666666666667 perc. fin"
## [1] "33.3333333333333 perc. fin"
## [1] "37.5 perc. fin"
## [1] "41.6666666666667 perc. fin"
## [1] "45.8333333333333 perc. fin"
## [1] "50 perc. fin"
## [1] "54.1666666666667 perc. fin"
## [1] "58.3333333333333 perc. fin"
## [1] "62.5 perc. fin"
## [1] "66.6666666666667 perc. fin"
## [1] "70.8333333333333 perc. fin"
## [1] "75 perc. fin"
## [1] "79.1666666666667 perc. fin"
## [1] "83.3333333333333 perc. fin"
## [1] "87.5 perc. fin"
## [1] "91.6666666666667 perc. fin"
## [1] "95.8333333333333 perc. fin"
## [1] "100 perc. fin"
## [1] "/Users/takanorihasegawa/Documents/GitHub/Neoantimon/vignettes/result.ID.SV2/sample_sv_bnd.txt.ID_SVFusion.HLACLASS2.4.peptide.txt"
## [1] "4.16666666666667 perc. fin"
## [1] "8.33333333333333 perc. fin"
## [1] "12.5 perc. fin"
## [1] "16.6666666666667 perc. fin"
## [1] "20.8333333333333 perc. fin"
## [1] "25 perc. fin"
## [1] "29.1666666666667 perc. fin"
## [1] "33.3333333333333 perc. fin"
## [1] "37.5 perc. fin"
## [1] "41.6666666666667 perc. fin"
## [1] "45.8333333333333 perc. fin"
## [1] "50 perc. fin"
## [1] "54.1666666666667 perc. fin"
## [1] "58.3333333333333 perc. fin"
## [1] "62.5 perc. fin"
## [1] "66.6666666666667 perc. fin"
## [1] "70.8333333333333 perc. fin"
## [1] "75 perc. fin"
## [1] "79.1666666666667 perc. fin"
## [1] "83.3333333333333 perc. fin"
## [1] "87.5 perc. fin"
## [1] "91.6666666666667 perc. fin"
## [1] "95.8333333333333 perc. fin"
## [1] "100 perc. fin"
## [1] "/Users/takanorihasegawa/Documents/GitHub/Neoantimon/vignettes/result.ID.SV2/sample_sv_bnd.txt.ID_SVFusion.HLACLASS2.5.peptide.txt"
## [1] "4.16666666666667 perc. fin"
## [1] "8.33333333333333 perc. fin"
## [1] "12.5 perc. fin"
## [1] "16.6666666666667 perc. fin"
## [1] "20.8333333333333 perc. fin"
## [1] "25 perc. fin"
## [1] "29.1666666666667 perc. fin"
## [1] "33.3333333333333 perc. fin"
## [1] "37.5 perc. fin"
## [1] "41.6666666666667 perc. fin"
## [1] "45.8333333333333 perc. fin"
## [1] "50 perc. fin"
## [1] "54.1666666666667 perc. fin"
## [1] "58.3333333333333 perc. fin"
## [1] "62.5 perc. fin"
## [1] "66.6666666666667 perc. fin"
## [1] "70.8333333333333 perc. fin"
## [1] "75 perc. fin"
## [1] "79.1666666666667 perc. fin"
## [1] "83.3333333333333 perc. fin"
## [1] "87.5 perc. fin"
## [1] "91.6666666666667 perc. fin"
## [1] "95.8333333333333 perc. fin"
## [1] "100 perc. fin"
## [1] "/Users/takanorihasegawa/Documents/GitHub/Neoantimon/vignettes/result.ID.SV2/sample_sv_bnd.txt.ID_SVFusion.HLACLASS2.6.peptide.txt"
## [1] "4.16666666666667 perc. fin"
## [1] "8.33333333333333 perc. fin"
## [1] "12.5 perc. fin"
## [1] "16.6666666666667 perc. fin"
## [1] "20.8333333333333 perc. fin"
## [1] "25 perc. fin"
## [1] "29.1666666666667 perc. fin"
## [1] "33.3333333333333 perc. fin"
## [1] "37.5 perc. fin"
## [1] "41.6666666666667 perc. fin"
## [1] "45.8333333333333 perc. fin"
## [1] "50 perc. fin"
## [1] "54.1666666666667 perc. fin"
## [1] "58.3333333333333 perc. fin"
## [1] "62.5 perc. fin"
## [1] "66.6666666666667 perc. fin"
## [1] "70.8333333333333 perc. fin"
## [1] "75 perc. fin"
## [1] "79.1666666666667 perc. fin"
## [1] "83.3333333333333 perc. fin"
## [1] "87.5 perc. fin"
## [1] "91.6666666666667 perc. fin"
## [1] "95.8333333333333 perc. fin"
## [1] "100 perc. fin"
## [1] "/Users/takanorihasegawa/Documents/GitHub/Neoantimon/vignettes/result.ID.SV2/sample_sv_bnd.txt.ID_SVFusion.HLACLASS2.7.peptide.txt"
## [1] "4.16666666666667 perc. fin"
## [1] "8.33333333333333 perc. fin"
## [1] "12.5 perc. fin"
## [1] "16.6666666666667 perc. fin"
## [1] "20.8333333333333 perc. fin"
## [1] "25 perc. fin"
## [1] "29.1666666666667 perc. fin"
## [1] "33.3333333333333 perc. fin"
## [1] "37.5 perc. fin"
## [1] "41.6666666666667 perc. fin"
## [1] "45.8333333333333 perc. fin"
## [1] "50 perc. fin"
## [1] "54.1666666666667 perc. fin"
## [1] "58.3333333333333 perc. fin"
## [1] "62.5 perc. fin"
## [1] "66.6666666666667 perc. fin"
## [1] "70.8333333333333 perc. fin"
## [1] "75 perc. fin"
## [1] "79.1666666666667 perc. fin"
## [1] "83.3333333333333 perc. fin"
## [1] "87.5 perc. fin"
## [1] "91.6666666666667 perc. fin"
## [1] "95.8333333333333 perc. fin"
## [1] "100 perc. fin"
## [1] "/Users/takanorihasegawa/Documents/GitHub/Neoantimon/vignettes/result.ID.SV2/sample_sv_bnd.txt.ID_SVFusion.HLACLASS2.8.peptide.txt"
## [1] "4.16666666666667 perc. fin"
## [1] "8.33333333333333 perc. fin"
## [1] "12.5 perc. fin"
## [1] "16.6666666666667 perc. fin"
## [1] "20.8333333333333 perc. fin"
## [1] "25 perc. fin"
## [1] "29.1666666666667 perc. fin"
## [1] "33.3333333333333 perc. fin"
## [1] "37.5 perc. fin"
## [1] "41.6666666666667 perc. fin"
## [1] "45.8333333333333 perc. fin"
## [1] "50 perc. fin"
## [1] "54.1666666666667 perc. fin"
## [1] "58.3333333333333 perc. fin"
## [1] "62.5 perc. fin"
## [1] "66.6666666666667 perc. fin"
## [1] "70.8333333333333 perc. fin"
## [1] "75 perc. fin"
## [1] "79.1666666666667 perc. fin"
## [1] "83.3333333333333 perc. fin"
## [1] "87.5 perc. fin"
## [1] "91.6666666666667 perc. fin"
## [1] "95.8333333333333 perc. fin"
## [1] "100 perc. fin"
## [1] "/Users/takanorihasegawa/Documents/GitHub/Neoantimon/vignettes/result.ID.SV2/sample_sv_bnd.txt.ID_SVFusion.HLACLASS2.9.peptide.txt"
## [1] "4.16666666666667 perc. fin"
## [1] "8.33333333333333 perc. fin"
## [1] "12.5 perc. fin"
## [1] "16.6666666666667 perc. fin"
## [1] "20.8333333333333 perc. fin"
## [1] "25 perc. fin"
## [1] "29.1666666666667 perc. fin"
## [1] "33.3333333333333 perc. fin"
## [1] "37.5 perc. fin"
## [1] "41.6666666666667 perc. fin"
## [1] "45.8333333333333 perc. fin"
## [1] "50 perc. fin"
## [1] "54.1666666666667 perc. fin"
## [1] "58.3333333333333 perc. fin"
## [1] "62.5 perc. fin"
## [1] "66.6666666666667 perc. fin"
## [1] "70.8333333333333 perc. fin"
## [1] "75 perc. fin"
## [1] "79.1666666666667 perc. fin"
## [1] "83.3333333333333 perc. fin"
## [1] "87.5 perc. fin"
## [1] "91.6666666666667 perc. fin"
## [1] "95.8333333333333 perc. fin"
## [1] "100 perc. fin"
## [1] "Successfully Finished."
```

``` r
  Result_HLA2_SV <- CalculatePriorityScores(result = Result_HLA2_SV, useRNAvaf = FALSE)
  print(head(Result_HLA2_SV))
```

```
##      HLA         Pos Gene     Evaluated_Mutant_Peptide_Core
## [1,] "DRB1_0101" "1" "0_AAR2" "YNSWEVGPK"                  
## [2,] "DRB1_0101" "2" "0_AAR2" "YNSWEVGPK"                  
## [3,] "DRB1_0101" "3" "0_AAR2" "YNSWEVGPK"                  
## [4,] "DRB1_0101" "4" "0_AAR2" "WEVGPKFRE"                  
## [5,] "DRB1_0101" "5" "0_AAR2" "WEVGPKFRE"                  
## [6,] "DRB1_0101" "6" "0_AAR2" "VGPKFREQL"                  
##      Evaluated_Mutant_Peptide Mut_EL     Mut_Rank Chr  NM_ID                
## [1,] "GIDYNSWEVGPKFRE"        "0.002971" "35.75"  "20" "NM_015511_NM_015906"
## [2,] "IDYNSWEVGPKFREQ"        "0.000624" "54.09"  "20" "NM_015511_NM_015906"
## [3,] "DYNSWEVGPKFREQL"        "0.000203" "68.45"  "20" "NM_015511_NM_015906"
## [4,] "YNSWEVGPKFREQLK"        "0.000105" "76.22"  "20" "NM_015511_NM_015906"
## [5,] "NSWEVGPKFREQLKL"        "0.000042" "85.71"  "20" "NM_015511_NM_015906"
## [6,] "SWEVGPKFREQLKLF"        "0.000009" "95.50"  "20" "NM_015511_NM_015906"
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
##      P_I                 P_R                     P                      
## [1,] "0.994835338984919" "5.13804056470373e-74"  "5.13805365806756e-74" 
## [2,] "0.994835358345505" "7.69118996528602e-114" "7.69120956488771e-114"
## [3,] "0.994835361818295" "5.05416874225877e-145" "5.05418162189048e-145"
## [4,] "0.994835362626685" "6.78124000735739e-162" "6.78125728811676e-162"
## [5,] "0.994835363146364" "1.67508152720577e-182" "1.67508579584702e-182"
## [6,] "0.994835363418576" "9.23253707018305e-204" "9.23256059762808e-204"
```


``` r
  print(Export_Summary_IndelSV(Result_HLA2_SV, Mut_Rank_th = 5))
```

```
##              Num_All_Alteration        Num_Evaluated_Alteration 
##                              16                              16 
## Num_Alteration_Generating_NeoAg                 Num_All_Peptide 
##                               8                             588 
##           Num_Evaluated_Peptide    Num_Peptide_Generating_NeoAg 
##                             588                              30
```


``` r
  print(Export_Summary_IndelSV_perFragments(Result_HLA2_SV, Mut_Rank_th = 5))
```

```
##                                                     GIDYNSWEVGPKFREQLKLFCETCDRLTCRDCQLLEHK-0_AAR2
## Num_Peptide_Per_Pep                                                                        14.000
## Num_Cond_Peptide_Per_Pep                                                                   14.000
## Num_Rest_Peptide_Per_Pep                                                                    0.000
## Num_Rest_Peptide_Per_Pep / Num_Cond_Peptide_Per_Pep                                         0.000
## -logP                                                                                       0.874
##                                                     KRGDPHELRNIFLQLSAVQEAQLKRLEVTRPRVLGSREQGQVPRMARQPPPPWVHAAFLLCLLSLGGAI-1_SLC25A12
## Num_Peptide_Per_Pep                                                                                                           44.000
## Num_Cond_Peptide_Per_Pep                                                                                                      44.000
## Num_Rest_Peptide_Per_Pep                                                                                                       3.000
## Num_Rest_Peptide_Per_Pep / Num_Cond_Peptide_Per_Pep                                                                            0.068
## -logP                                                                                                                          1.588
##                                                     RKTMAEKRQLFIEMLYRLX-2_DTNB
## Num_Peptide_Per_Pep                                                      4.000
## Num_Cond_Peptide_Per_Pep                                                 4.000
## Num_Rest_Peptide_Per_Pep                                                 0.000
## Num_Rest_Peptide_Per_Pep / Num_Cond_Peptide_Per_Pep                      0.000
## -logP                                                                    0.437
##                                                     GHTDWLSDCCFHPSERCRCTLYGHTDSVNSIEFFPFS-3_SPAG16
## Num_Peptide_Per_Pep                                                                         13.000
## Num_Cond_Peptide_Per_Pep                                                                    13.000
## Num_Rest_Peptide_Per_Pep                                                                     0.000
## Num_Rest_Peptide_Per_Pep / Num_Cond_Peptide_Per_Pep                                          0.000
## -logP                                                                                        0.851
##                                                     MKSCKRNLSTTTGSSERX-4_TMCC1
## Num_Peptide_Per_Pep                                                      3.000
## Num_Cond_Peptide_Per_Pep                                                 3.000
## Num_Rest_Peptide_Per_Pep                                                 0.000
## Num_Rest_Peptide_Per_Pep / Num_Cond_Peptide_Per_Pep                      0.000
## -logP                                                                    0.414
##                                                     AAACLPGEEGTAERSCKRNLSTTTGSSERX-5_TMCC1
## Num_Peptide_Per_Pep                                                                  15.00
## Num_Cond_Peptide_Per_Pep                                                             15.00
## Num_Rest_Peptide_Per_Pep                                                              0.00
## Num_Rest_Peptide_Per_Pep / Num_Cond_Peptide_Per_Pep                                   0.00
## -logP                                                                                 0.69
##                                                     MVQRFSLRRQLSKSCKRNLSTTTGSSERX-6_TMCC1
## Num_Peptide_Per_Pep                                                                14.000
## Num_Cond_Peptide_Per_Pep                                                           14.000
## Num_Rest_Peptide_Per_Pep                                                            0.000
## Num_Rest_Peptide_Per_Pep / Num_Cond_Peptide_Per_Pep                                 0.000
## -logP                                                                               0.667
##                                                     PMQRKGSYCWDFEESCKRNLSTTTGSSERX-7_TMCC1
## Num_Peptide_Per_Pep                                                                  15.00
## Num_Cond_Peptide_Per_Pep                                                             15.00
## Num_Rest_Peptide_Per_Pep                                                              0.00
## Num_Rest_Peptide_Per_Pep / Num_Cond_Peptide_Per_Pep                                   0.00
## -logP                                                                                 0.69
##                                                     HNIRPKPFVIPGRSRTLECKQQKIILLRREMKETX-8_PDLIM3
## Num_Peptide_Per_Pep                                                                       20.000
## Num_Cond_Peptide_Per_Pep                                                                  20.000
## Num_Rest_Peptide_Per_Pep                                                                   2.000
## Num_Rest_Peptide_Per_Pep / Num_Cond_Peptide_Per_Pep                                        0.100
## -logP                                                                                      0.805
##                                                     LHGQLRGLIPSSPQKTLECKQQKIILLRREMKETX-9_PDLIM3
## Num_Peptide_Per_Pep                                                                       20.000
## Num_Cond_Peptide_Per_Pep                                                                  20.000
## Num_Rest_Peptide_Per_Pep                                                                   2.000
## Num_Rest_Peptide_Per_Pep / Num_Cond_Peptide_Per_Pep                                        0.100
## -logP                                                                                      0.805
##                                                     IKGDQYIVDMANTKRTLLWTRPMTVILGKPFCVMPKQQKTAHIGFLQHIPRLSPKPCLPKLNLMMRKQRMSQNGKNVKFEESHLRAVCMSGRGMGQVWVFFLCSX-10_WDR70
## Num_Peptide_Per_Pep                                                                                                                                             90.000
## Num_Cond_Peptide_Per_Pep                                                                                                                                        90.000
## Num_Rest_Peptide_Per_Pep                                                                                                                                        16.000
## Num_Rest_Peptide_Per_Pep / Num_Cond_Peptide_Per_Pep                                                                                                              0.178
## -logP                                                                                                                                                            2.416
##                                                     AALCPASRALEEKKGHELLGGPSLGAPRPGSQERTGENTAACQDHRFWAGQTAGCGRERIPCRRRQSAYQVDGIGINFTQNLYPPEX-11_EGFR
## Num_Peptide_Per_Pep                                                                                                                          72.000
## Num_Cond_Peptide_Per_Pep                                                                                                                     72.000
## Num_Rest_Peptide_Per_Pep                                                                                                                      1.000
## Num_Rest_Peptide_Per_Pep / Num_Cond_Peptide_Per_Pep                                                                                           0.014
## -logP                                                                                                                                         2.002
##                                                     AILTPMPEGLSQQQDMLKNTSKGHPDRLPLQMALTEL-12_ARHGEF10
## Num_Peptide_Per_Pep                                                                            13.000
## Num_Cond_Peptide_Per_Pep                                                                       13.000
## Num_Rest_Peptide_Per_Pep                                                                        1.000
## Num_Rest_Peptide_Per_Pep / Num_Cond_Peptide_Per_Pep                                             0.077
## -logP                                                                                           0.851
##                                                     LKPGVQGKKVRKPSPTTSASSSSWTSAGTSWPAGRRTGTATSGTATTTSVWPGCGTRMWSTQWSSVPRSRSCCSRPATTPPSKPGAPHAPCASSRHLAHGLAPSSPGLPARGAEVCWVHWSHRDPLRTSPGSVAFSRAGEVEMLIAVTPX-13_NOTCH1
## Num_Peptide_Per_Pep                                                                                                                                                                                          135.000
## Num_Cond_Peptide_Per_Pep                                                                                                                                                                                     135.000
## Num_Rest_Peptide_Per_Pep                                                                                                                                                                                       1.000
## Num_Rest_Peptide_Per_Pep / Num_Cond_Peptide_Per_Pep                                                                                                                                                            0.007
## -logP                                                                                                                                                                                                          3.451
##                                                     ESSRGRGNNNRKGREVTPL-14_UBAP2
## Num_Peptide_Per_Pep                                                        5.000
## Num_Cond_Peptide_Per_Pep                                                   5.000
## Num_Rest_Peptide_Per_Pep                                                   0.000
## Num_Rest_Peptide_Per_Pep / Num_Cond_Peptide_Per_Pep                        0.000
## -logP                                                                      0.437
##                                                     WYIPSSQLAAKLLHMLTMRMLSKSATGRWX-15_RASGRP2
## Num_Peptide_Per_Pep                                                                     15.00
## Num_Cond_Peptide_Per_Pep                                                                15.00
## Num_Rest_Peptide_Per_Pep                                                                 3.00
## Num_Rest_Peptide_Per_Pep / Num_Cond_Peptide_Per_Pep                                      0.20
## -logP                                                                                    0.69
##                                                     WLRKSPPEKKLRRYPPGQGATIGVDFMIKTVEINGEKV-16_GAB2
## Num_Peptide_Per_Pep                                                                         14.000
## Num_Cond_Peptide_Per_Pep                                                                    14.000
## Num_Rest_Peptide_Per_Pep                                                                     0.000
## Num_Rest_Peptide_Per_Pep / Num_Cond_Peptide_Per_Pep                                          0.000
## -logP                                                                                        0.874
##                                                     GTMSVIVLPQPDEVLNLVQSYVTLRVPLYVSYVFHSP-17_MTRF1
## Num_Peptide_Per_Pep                                                                         13.000
## Num_Cond_Peptide_Per_Pep                                                                    13.000
## Num_Rest_Peptide_Per_Pep                                                                     0.000
## Num_Rest_Peptide_Per_Pep / Num_Cond_Peptide_Per_Pep                                          0.000
## -logP                                                                                        0.851
##                                                     IVLWKCLRKPVLNDRNLQLHTDKGSFLKEKNKKLKKKNK-18_EVI2B
## Num_Peptide_Per_Pep                                                                           15.000
## Num_Cond_Peptide_Per_Pep                                                                      15.000
## Num_Rest_Peptide_Per_Pep                                                                       0.000
## Num_Rest_Peptide_Per_Pep / Num_Cond_Peptide_Per_Pep                                            0.000
## -logP                                                                                          0.897
##                                                     RCRTFLQSTDPGDTGGVFEGSTLCPSDPAFSRTDDPX-19_IKZF3
## Num_Peptide_Per_Pep                                                                         22.000
## Num_Cond_Peptide_Per_Pep                                                                    22.000
## Num_Rest_Peptide_Per_Pep                                                                     0.000
## Num_Rest_Peptide_Per_Pep / Num_Cond_Peptide_Per_Pep                                          0.000
## -logP                                                                                        0.851
##                                                     QRRDALTGHLRTHSGGVFEGSTLCPSDPAFSRTDDPX-20_IKZF3
## Num_Peptide_Per_Pep                                                                         22.000
## Num_Cond_Peptide_Per_Pep                                                                    22.000
## Num_Rest_Peptide_Per_Pep                                                                     0.000
## Num_Rest_Peptide_Per_Pep / Num_Cond_Peptide_Per_Pep                                          0.000
## -logP                                                                                        0.851
##                                                     DSGEGPANEDEDIGGGVFEGSTLCPSDPAFSRTDDPX-21_IKZF3
## Num_Peptide_Per_Pep                                                                         22.000
## Num_Cond_Peptide_Per_Pep                                                                    22.000
## Num_Rest_Peptide_Per_Pep                                                                     0.000
## Num_Rest_Peptide_Per_Pep / Num_Cond_Peptide_Per_Pep                                          0.000
## -logP                                                                                        0.851
##                                                     SFNVLMVHKRSHTGGVFEGSTLCPSDPAFSRTDDPX-22_IKZF3
## Num_Peptide_Per_Pep                                                                        21.000
## Num_Cond_Peptide_Per_Pep                                                                   21.000
## Num_Rest_Peptide_Per_Pep                                                                    2.000
## Num_Rest_Peptide_Per_Pep / Num_Cond_Peptide_Per_Pep                                         0.095
## -logP                                                                                       0.828
##                                                     ISFNVLMVHKRSHTGGVFEGSTLCPSDPAFSRTDDPX-23_IKZF3
## Num_Peptide_Per_Pep                                                                         22.000
## Num_Cond_Peptide_Per_Pep                                                                    22.000
## Num_Rest_Peptide_Per_Pep                                                                     3.000
## Num_Rest_Peptide_Per_Pep / Num_Cond_Peptide_Per_Pep                                          0.136
## -logP                                                                                        0.851
```

To calculate the binding affinity of neoantigen candaites, which are directly generated from RNA sequences, to HLA ClassI. The peptides included in the original genes ("NM_003998", "NM_001165412") are removed from the results.


``` r
   Result_HLA1_Seq <- MainSeqFragmentClass1(input_sequence = "atggcagaagatgatccatatttgggaaggcctgaaaaaatgtttcatttggatccttctttgactcatacaatatttaatc",
                                           file_name_in_hla_table = "sample",
                                           hla_file = "./../data/sample_hla_table_c1.txt",
                                           hmdir = getwd(),
                                           job_id = "NO_job_id",
                                           refflat_file  = "./../lib/refFlat.grch37.txt",
                                           refmrna_file = "./../lib/refMrna.grch37.fa",
                                           netMHCpan_dir = "./../lib/netMHCpan-4.1/netMHCpan",
                                           reference_nm_id = c("NM_003998", "NM_001165412"))
```

```
## [1] "Sequence: atggcagaagatgatccatatttgggaaggcctgaaaaaatgtttcatttggatccttctttgactcatacaatatttaatc"
## [1] "HLAtype: A*02:01" "HLAtype: A*32:01" "HLAtype: B*15:17" "HLAtype: B*51:01"
## [5] "HLAtype: C*07:01" "HLAtype: C*15:02"
## [1] "refFLAT: ./../lib/refFlat.grch37.txt"
## [1] "refMrna: ./../lib/refMrna.grch37.fa"
## [1] "Wt-NM_ID:  NM_003998"    "Wt-NM_ID:  NM_001165412"
## [1] "Wt-Gene Symbol: NA"
## [1] "Start Position of Reading Frame is: 1"
## [1] "Peptide Successfully Generated!!"
## [1] "Executing netMHCpan to result.NO_job_id.SeqFragment1"
## [1] "./../lib/netMHCpan-4.1/netMHCpan -BA  -l 8,9,10,11,12,13 -f result.NO_job_id.SeqFragment1/NO_job_id_SeqFragment.peptide.fasta -a HLA-A02:01 > result.NO_job_id.SeqFragment1/NO_job_id_SeqFragment.HLACLASS1.1.peptide.txt"
## [1] "./../lib/netMHCpan-4.1/netMHCpan -BA  -l 8,9,10,11,12,13 -f result.NO_job_id.SeqFragment1/NO_job_id_SeqFragment.peptide.fasta -a HLA-A32:01 > result.NO_job_id.SeqFragment1/NO_job_id_SeqFragment.HLACLASS1.2.peptide.txt"
## [1] "./../lib/netMHCpan-4.1/netMHCpan -BA  -l 8,9,10,11,12,13 -f result.NO_job_id.SeqFragment1/NO_job_id_SeqFragment.peptide.fasta -a HLA-B15:17 > result.NO_job_id.SeqFragment1/NO_job_id_SeqFragment.HLACLASS1.3.peptide.txt"
## [1] "./../lib/netMHCpan-4.1/netMHCpan -BA  -l 8,9,10,11,12,13 -f result.NO_job_id.SeqFragment1/NO_job_id_SeqFragment.peptide.fasta -a HLA-B51:01 > result.NO_job_id.SeqFragment1/NO_job_id_SeqFragment.HLACLASS1.4.peptide.txt"
## [1] "./../lib/netMHCpan-4.1/netMHCpan -BA  -l 8,9,10,11,12,13 -f result.NO_job_id.SeqFragment1/NO_job_id_SeqFragment.peptide.fasta -a HLA-C07:01 > result.NO_job_id.SeqFragment1/NO_job_id_SeqFragment.HLACLASS1.5.peptide.txt"
## [1] "./../lib/netMHCpan-4.1/netMHCpan -BA  -l 8,9,10,11,12,13 -f result.NO_job_id.SeqFragment1/NO_job_id_SeqFragment.peptide.fasta -a HLA-C15:02 > result.NO_job_id.SeqFragment1/NO_job_id_SeqFragment.HLACLASS1.6.peptide.txt"
## [1] "Merging Results..."
## [1] "/Users/takanorihasegawa/Documents/GitHub/Neoantimon/vignettes/result.NO_job_id.SeqFragment1/NO_job_id_SeqFragment.HLACLASS1.1.peptide.txt"
## [1] "100 perc. fin"
## [1] "/Users/takanorihasegawa/Documents/GitHub/Neoantimon/vignettes/result.NO_job_id.SeqFragment1/NO_job_id_SeqFragment.HLACLASS1.2.peptide.txt"
## [1] "100 perc. fin"
## [1] "/Users/takanorihasegawa/Documents/GitHub/Neoantimon/vignettes/result.NO_job_id.SeqFragment1/NO_job_id_SeqFragment.HLACLASS1.3.peptide.txt"
## [1] "100 perc. fin"
## [1] "/Users/takanorihasegawa/Documents/GitHub/Neoantimon/vignettes/result.NO_job_id.SeqFragment1/NO_job_id_SeqFragment.HLACLASS1.4.peptide.txt"
## [1] "100 perc. fin"
## [1] "/Users/takanorihasegawa/Documents/GitHub/Neoantimon/vignettes/result.NO_job_id.SeqFragment1/NO_job_id_SeqFragment.HLACLASS1.5.peptide.txt"
## [1] "100 perc. fin"
## [1] "/Users/takanorihasegawa/Documents/GitHub/Neoantimon/vignettes/result.NO_job_id.SeqFragment1/NO_job_id_SeqFragment.HLACLASS1.6.peptide.txt"
## [1] "100 perc. fin"
## [1] "Successfully Finished."
```

``` r
  Result_HLA1_Seq <- CalculatePriorityScores(result = Result_HLA1_Seq, useRNAvaf = FALSE)
```

```
## Warning in CalculatePriorityScores(result = Result_HLA1_Seq, useRNAvaf =
## FALSE): NAs introduced by coercion
```

``` r
  print(head(Result_HLA1_Seq))
```

```
##      HLA           Pos Gene           Evaluated_Mutant_Peptide_Core
## [1,] "HLA-A*02:01" "1" "0_atggcagaag" "MAEDDPYLGRPEK"              
## [2,] "HLA-A*02:01" "2" "0_atggcagaag" "AEDDPYLGRPEKM"              
## [3,] "HLA-A*02:01" "3" "0_atggcagaag" "EDDPYLGRPEKMF"              
## [4,] "HLA-A*02:01" "4" "0_atggcagaag" "DDPYLGRPEKMFH"              
## [5,] "HLA-A*02:01" "5" "0_atggcagaag" "YLGRPEKMFHL"                
## [6,] "HLA-A*02:01" "6" "0_atggcagaag" "YLGRPEKMFHL"                
##      Evaluated_Mutant_Peptide Mut_EL      Mut_Rank Chr NM_ID ReadingFrame
## [1,] "MAEDDPYLGRPEK"          "0.0000280" "59.545" "0" ""    "1"         
## [2,] "AEDDPYLGRPEKM"          "0.0001790" "33.938" "0" ""    "1"         
## [3,] "EDDPYLGRPEKMF"          "0.0000050" "82.500" "0" ""    "1"         
## [4,] "DDPYLGRPEKMFH"          "0.0000010" "95.000" "0" ""    "1"         
## [5,] "DPYLGRPEKMFHL"          "0.0001750" "34.214" "0" ""    "1"         
## [6,] "PYLGRPEKMFHLD"          "0.0000040" "85.000" "0" ""    "1"         
##      SequenceNumber Chrs        NM_IDs                   GeneIDs      
## [1,] "1"            "chr4;chr4" "NM_003998;NM_001165412" "NFKB1;NFKB1"
## [2,] "1"            "chr4;chr4" "NM_003998;NM_001165412" "NFKB1;NFKB1"
## [3,] "1"            "chr4;chr4" "NM_003998;NM_001165412" "NFKB1;NFKB1"
## [4,] "1"            "chr4;chr4" "NM_003998;NM_001165412" "NFKB1;NFKB1"
## [5,] "1"            "chr4;chr4" "NM_003998;NM_001165412" "NFKB1;NFKB1"
## [6,] "1"            "chr4;chr4" "NM_003998;NM_001165412" "NFKB1;NFKB1"
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
## [1,] "MAEDDPYLGRPEQMFHLDPSLTHTIFNPEVFQPQMALPTADGPYLQILEQPKQRGFRFRYVCEGPSHGGLPGASSEKNKKSYPQVKICNYVGPAKVIVQLVTNGKNIHLHAHSLVGKHCEDGICTVTAGPKDMVVGFANLGILHVTKKKVFETLEARMTEACIRGYNPGLLVHPDLAYLQAEGGGDRQLGDREKELIRQAALQQTKEMDLSVVRLMFTAFLPDSTGSFTRRLEPVVSDAIYDSKAPNASNLKIVRMDRTAGCVTGGEEIYLLCDKVQKDDIQIRFYEEEENGGVWEGFGDFSPTDVHRQFAIVFKTPKYKDINITKPASVFVQLRRKSDLETSEPKPFLYYPEIKDKEEVQRKRQKLMPNFSDSFGGGSGAGAGGGGMFGSGGGGGGTGSTGPGYSFPHYGFPTYGGITFHPGTTKSNAGMKHGTMDTESKKDPEGCDKSDDKNTVNLFGKVIETTEQDQEPSEATVGNGEVTLTYATGTKEESAGVQDNLFLEKAMQLAKRHANALFDYAVTGDVKMLLAVQRHLTAVQDENGDSVLHLAIIHLHSQLVRDLLEVTSGLISDDIINMRNDLYQTPLHLAVITKQEDVVEDLLRAGADLSLLDRLGNSVLHLAAKEGHDKVLSILLKHKKAALLLDHPNGDGLNAIHLAMMSNSLPCLLLLVAAGADVNAQEQKSGRTALHLAVEHDNISLAGCLLLEGDAHVDSTTYDGTTPLHIAAGRGSTRLAALLKAAGADPLVENFEPLYDLDDSWENAGEDEGVVPGTTPLDMATSWQVFDILNGKPYEPEFTSDDLLAQGDMKQLAEDVKLQLYKLLEIPDPDKNWATLAQKLGLGILNNAFRLSPAPSKTLMDNYEVSGGTVRELVEALRQMGYTEAIEVIQAASSPVKTTSQAHSLPLSPASTRQQIDELRDSDSVCDSGVETSFRKLSFTESLTSGASLLTLNKMPHDYGQEGPLEGKIXXMAEDDPYLGRPEQMFHLDPSLTHTIFNPEVFQPQMALPTDGPYLQILEQPKQRGFRFRYVCEGPSHGGLPGASSEKNKKSYPQVKICNYVGPAKVIVQLVTNGKNIHLHAHSLVGKHCEDGICTVTAGPKDMVVGFANLGILHVTKKKVFETLEARMTEACIRGYNPGLLVHPDLAYLQAEGGGDRQLGDREKELIRQAALQQTKEMDLSVVRLMFTAFLPDSTGSFTRRLEPVVSDAIYDSKAPNASNLKIVRMDRTAGCVTGGEEIYLLCDKVQKDDIQIRFYEEEENGGVWEGFGDFSPTDVHRQFAIVFKTPKYKDINITKPASVFVQLRRKSDLETSEPKPFLYYPEIKDKEEVQRKRQKLMPNFSDSFGGGSGAGAGGGGMFGSGGGGGGTGSTGPGYSFPHYGFPTYGGITFHPGTTKSNAGMKHGTMDTESKKDPEGCDKSDDKNTVNLFGKVIETTEQDQEPSEATVGNGEVTLTYATGTKEESAGVQDNLFLEKAMQLAKRHANALFDYAVTGDVKMLLAVQRHLTAVQDENGDSVLHLAIIHLHSQLVRDLLEVTSGLISDDIINMRNDLYQTPLHLAVITKQEDVVEDLLRAGADLSLLDRLGNSVLHLAAKEGHDKVLSILLKHKKAALLLDHPNGDGLNAIHLAMMSNSLPCLLLLVAAGADVNAQEQKSGRTALHLAVEHDNISLAGCLLLEGDAHVDSTTYDGTTPLHIAAGRGSTRLAALLKAAGADPLVENFEPLYDLDDSWENAGEDEGVVPGTTPLDMATSWQVFDILNGKPYEPEFTSDDLLAQGDMKQLAEDVKLQLYKLLEIPDPDKNWATLAQKLGLGILNNAFRLSPAPSKTLMDNYEVSGGTVRELVEALRQMGYTEAIEVIQAASSPVKTTSQAHSLPLSPASTRQQIDELRDSDSVCDSGVETSFRKLSFTESLTSGASLLTLNKMPHDYGQEGPLEGKIX"
## [2,] "MAEDDPYLGRPEQMFHLDPSLTHTIFNPEVFQPQMALPTADGPYLQILEQPKQRGFRFRYVCEGPSHGGLPGASSEKNKKSYPQVKICNYVGPAKVIVQLVTNGKNIHLHAHSLVGKHCEDGICTVTAGPKDMVVGFANLGILHVTKKKVFETLEARMTEACIRGYNPGLLVHPDLAYLQAEGGGDRQLGDREKELIRQAALQQTKEMDLSVVRLMFTAFLPDSTGSFTRRLEPVVSDAIYDSKAPNASNLKIVRMDRTAGCVTGGEEIYLLCDKVQKDDIQIRFYEEEENGGVWEGFGDFSPTDVHRQFAIVFKTPKYKDINITKPASVFVQLRRKSDLETSEPKPFLYYPEIKDKEEVQRKRQKLMPNFSDSFGGGSGAGAGGGGMFGSGGGGGGTGSTGPGYSFPHYGFPTYGGITFHPGTTKSNAGMKHGTMDTESKKDPEGCDKSDDKNTVNLFGKVIETTEQDQEPSEATVGNGEVTLTYATGTKEESAGVQDNLFLEKAMQLAKRHANALFDYAVTGDVKMLLAVQRHLTAVQDENGDSVLHLAIIHLHSQLVRDLLEVTSGLISDDIINMRNDLYQTPLHLAVITKQEDVVEDLLRAGADLSLLDRLGNSVLHLAAKEGHDKVLSILLKHKKAALLLDHPNGDGLNAIHLAMMSNSLPCLLLLVAAGADVNAQEQKSGRTALHLAVEHDNISLAGCLLLEGDAHVDSTTYDGTTPLHIAAGRGSTRLAALLKAAGADPLVENFEPLYDLDDSWENAGEDEGVVPGTTPLDMATSWQVFDILNGKPYEPEFTSDDLLAQGDMKQLAEDVKLQLYKLLEIPDPDKNWATLAQKLGLGILNNAFRLSPAPSKTLMDNYEVSGGTVRELVEALRQMGYTEAIEVIQAASSPVKTTSQAHSLPLSPASTRQQIDELRDSDSVCDSGVETSFRKLSFTESLTSGASLLTLNKMPHDYGQEGPLEGKIXXMAEDDPYLGRPEQMFHLDPSLTHTIFNPEVFQPQMALPTDGPYLQILEQPKQRGFRFRYVCEGPSHGGLPGASSEKNKKSYPQVKICNYVGPAKVIVQLVTNGKNIHLHAHSLVGKHCEDGICTVTAGPKDMVVGFANLGILHVTKKKVFETLEARMTEACIRGYNPGLLVHPDLAYLQAEGGGDRQLGDREKELIRQAALQQTKEMDLSVVRLMFTAFLPDSTGSFTRRLEPVVSDAIYDSKAPNASNLKIVRMDRTAGCVTGGEEIYLLCDKVQKDDIQIRFYEEEENGGVWEGFGDFSPTDVHRQFAIVFKTPKYKDINITKPASVFVQLRRKSDLETSEPKPFLYYPEIKDKEEVQRKRQKLMPNFSDSFGGGSGAGAGGGGMFGSGGGGGGTGSTGPGYSFPHYGFPTYGGITFHPGTTKSNAGMKHGTMDTESKKDPEGCDKSDDKNTVNLFGKVIETTEQDQEPSEATVGNGEVTLTYATGTKEESAGVQDNLFLEKAMQLAKRHANALFDYAVTGDVKMLLAVQRHLTAVQDENGDSVLHLAIIHLHSQLVRDLLEVTSGLISDDIINMRNDLYQTPLHLAVITKQEDVVEDLLRAGADLSLLDRLGNSVLHLAAKEGHDKVLSILLKHKKAALLLDHPNGDGLNAIHLAMMSNSLPCLLLLVAAGADVNAQEQKSGRTALHLAVEHDNISLAGCLLLEGDAHVDSTTYDGTTPLHIAAGRGSTRLAALLKAAGADPLVENFEPLYDLDDSWENAGEDEGVVPGTTPLDMATSWQVFDILNGKPYEPEFTSDDLLAQGDMKQLAEDVKLQLYKLLEIPDPDKNWATLAQKLGLGILNNAFRLSPAPSKTLMDNYEVSGGTVRELVEALRQMGYTEAIEVIQAASSPVKTTSQAHSLPLSPASTRQQIDELRDSDSVCDSGVETSFRKLSFTESLTSGASLLTLNKMPHDYGQEGPLEGKIX"
## [3,] "MAEDDPYLGRPEQMFHLDPSLTHTIFNPEVFQPQMALPTADGPYLQILEQPKQRGFRFRYVCEGPSHGGLPGASSEKNKKSYPQVKICNYVGPAKVIVQLVTNGKNIHLHAHSLVGKHCEDGICTVTAGPKDMVVGFANLGILHVTKKKVFETLEARMTEACIRGYNPGLLVHPDLAYLQAEGGGDRQLGDREKELIRQAALQQTKEMDLSVVRLMFTAFLPDSTGSFTRRLEPVVSDAIYDSKAPNASNLKIVRMDRTAGCVTGGEEIYLLCDKVQKDDIQIRFYEEEENGGVWEGFGDFSPTDVHRQFAIVFKTPKYKDINITKPASVFVQLRRKSDLETSEPKPFLYYPEIKDKEEVQRKRQKLMPNFSDSFGGGSGAGAGGGGMFGSGGGGGGTGSTGPGYSFPHYGFPTYGGITFHPGTTKSNAGMKHGTMDTESKKDPEGCDKSDDKNTVNLFGKVIETTEQDQEPSEATVGNGEVTLTYATGTKEESAGVQDNLFLEKAMQLAKRHANALFDYAVTGDVKMLLAVQRHLTAVQDENGDSVLHLAIIHLHSQLVRDLLEVTSGLISDDIINMRNDLYQTPLHLAVITKQEDVVEDLLRAGADLSLLDRLGNSVLHLAAKEGHDKVLSILLKHKKAALLLDHPNGDGLNAIHLAMMSNSLPCLLLLVAAGADVNAQEQKSGRTALHLAVEHDNISLAGCLLLEGDAHVDSTTYDGTTPLHIAAGRGSTRLAALLKAAGADPLVENFEPLYDLDDSWENAGEDEGVVPGTTPLDMATSWQVFDILNGKPYEPEFTSDDLLAQGDMKQLAEDVKLQLYKLLEIPDPDKNWATLAQKLGLGILNNAFRLSPAPSKTLMDNYEVSGGTVRELVEALRQMGYTEAIEVIQAASSPVKTTSQAHSLPLSPASTRQQIDELRDSDSVCDSGVETSFRKLSFTESLTSGASLLTLNKMPHDYGQEGPLEGKIXXMAEDDPYLGRPEQMFHLDPSLTHTIFNPEVFQPQMALPTDGPYLQILEQPKQRGFRFRYVCEGPSHGGLPGASSEKNKKSYPQVKICNYVGPAKVIVQLVTNGKNIHLHAHSLVGKHCEDGICTVTAGPKDMVVGFANLGILHVTKKKVFETLEARMTEACIRGYNPGLLVHPDLAYLQAEGGGDRQLGDREKELIRQAALQQTKEMDLSVVRLMFTAFLPDSTGSFTRRLEPVVSDAIYDSKAPNASNLKIVRMDRTAGCVTGGEEIYLLCDKVQKDDIQIRFYEEEENGGVWEGFGDFSPTDVHRQFAIVFKTPKYKDINITKPASVFVQLRRKSDLETSEPKPFLYYPEIKDKEEVQRKRQKLMPNFSDSFGGGSGAGAGGGGMFGSGGGGGGTGSTGPGYSFPHYGFPTYGGITFHPGTTKSNAGMKHGTMDTESKKDPEGCDKSDDKNTVNLFGKVIETTEQDQEPSEATVGNGEVTLTYATGTKEESAGVQDNLFLEKAMQLAKRHANALFDYAVTGDVKMLLAVQRHLTAVQDENGDSVLHLAIIHLHSQLVRDLLEVTSGLISDDIINMRNDLYQTPLHLAVITKQEDVVEDLLRAGADLSLLDRLGNSVLHLAAKEGHDKVLSILLKHKKAALLLDHPNGDGLNAIHLAMMSNSLPCLLLLVAAGADVNAQEQKSGRTALHLAVEHDNISLAGCLLLEGDAHVDSTTYDGTTPLHIAAGRGSTRLAALLKAAGADPLVENFEPLYDLDDSWENAGEDEGVVPGTTPLDMATSWQVFDILNGKPYEPEFTSDDLLAQGDMKQLAEDVKLQLYKLLEIPDPDKNWATLAQKLGLGILNNAFRLSPAPSKTLMDNYEVSGGTVRELVEALRQMGYTEAIEVIQAASSPVKTTSQAHSLPLSPASTRQQIDELRDSDSVCDSGVETSFRKLSFTESLTSGASLLTLNKMPHDYGQEGPLEGKIX"
## [4,] "MAEDDPYLGRPEQMFHLDPSLTHTIFNPEVFQPQMALPTADGPYLQILEQPKQRGFRFRYVCEGPSHGGLPGASSEKNKKSYPQVKICNYVGPAKVIVQLVTNGKNIHLHAHSLVGKHCEDGICTVTAGPKDMVVGFANLGILHVTKKKVFETLEARMTEACIRGYNPGLLVHPDLAYLQAEGGGDRQLGDREKELIRQAALQQTKEMDLSVVRLMFTAFLPDSTGSFTRRLEPVVSDAIYDSKAPNASNLKIVRMDRTAGCVTGGEEIYLLCDKVQKDDIQIRFYEEEENGGVWEGFGDFSPTDVHRQFAIVFKTPKYKDINITKPASVFVQLRRKSDLETSEPKPFLYYPEIKDKEEVQRKRQKLMPNFSDSFGGGSGAGAGGGGMFGSGGGGGGTGSTGPGYSFPHYGFPTYGGITFHPGTTKSNAGMKHGTMDTESKKDPEGCDKSDDKNTVNLFGKVIETTEQDQEPSEATVGNGEVTLTYATGTKEESAGVQDNLFLEKAMQLAKRHANALFDYAVTGDVKMLLAVQRHLTAVQDENGDSVLHLAIIHLHSQLVRDLLEVTSGLISDDIINMRNDLYQTPLHLAVITKQEDVVEDLLRAGADLSLLDRLGNSVLHLAAKEGHDKVLSILLKHKKAALLLDHPNGDGLNAIHLAMMSNSLPCLLLLVAAGADVNAQEQKSGRTALHLAVEHDNISLAGCLLLEGDAHVDSTTYDGTTPLHIAAGRGSTRLAALLKAAGADPLVENFEPLYDLDDSWENAGEDEGVVPGTTPLDMATSWQVFDILNGKPYEPEFTSDDLLAQGDMKQLAEDVKLQLYKLLEIPDPDKNWATLAQKLGLGILNNAFRLSPAPSKTLMDNYEVSGGTVRELVEALRQMGYTEAIEVIQAASSPVKTTSQAHSLPLSPASTRQQIDELRDSDSVCDSGVETSFRKLSFTESLTSGASLLTLNKMPHDYGQEGPLEGKIXXMAEDDPYLGRPEQMFHLDPSLTHTIFNPEVFQPQMALPTDGPYLQILEQPKQRGFRFRYVCEGPSHGGLPGASSEKNKKSYPQVKICNYVGPAKVIVQLVTNGKNIHLHAHSLVGKHCEDGICTVTAGPKDMVVGFANLGILHVTKKKVFETLEARMTEACIRGYNPGLLVHPDLAYLQAEGGGDRQLGDREKELIRQAALQQTKEMDLSVVRLMFTAFLPDSTGSFTRRLEPVVSDAIYDSKAPNASNLKIVRMDRTAGCVTGGEEIYLLCDKVQKDDIQIRFYEEEENGGVWEGFGDFSPTDVHRQFAIVFKTPKYKDINITKPASVFVQLRRKSDLETSEPKPFLYYPEIKDKEEVQRKRQKLMPNFSDSFGGGSGAGAGGGGMFGSGGGGGGTGSTGPGYSFPHYGFPTYGGITFHPGTTKSNAGMKHGTMDTESKKDPEGCDKSDDKNTVNLFGKVIETTEQDQEPSEATVGNGEVTLTYATGTKEESAGVQDNLFLEKAMQLAKRHANALFDYAVTGDVKMLLAVQRHLTAVQDENGDSVLHLAIIHLHSQLVRDLLEVTSGLISDDIINMRNDLYQTPLHLAVITKQEDVVEDLLRAGADLSLLDRLGNSVLHLAAKEGHDKVLSILLKHKKAALLLDHPNGDGLNAIHLAMMSNSLPCLLLLVAAGADVNAQEQKSGRTALHLAVEHDNISLAGCLLLEGDAHVDSTTYDGTTPLHIAAGRGSTRLAALLKAAGADPLVENFEPLYDLDDSWENAGEDEGVVPGTTPLDMATSWQVFDILNGKPYEPEFTSDDLLAQGDMKQLAEDVKLQLYKLLEIPDPDKNWATLAQKLGLGILNNAFRLSPAPSKTLMDNYEVSGGTVRELVEALRQMGYTEAIEVIQAASSPVKTTSQAHSLPLSPASTRQQIDELRDSDSVCDSGVETSFRKLSFTESLTSGASLLTLNKMPHDYGQEGPLEGKIX"
## [5,] "MAEDDPYLGRPEQMFHLDPSLTHTIFNPEVFQPQMALPTADGPYLQILEQPKQRGFRFRYVCEGPSHGGLPGASSEKNKKSYPQVKICNYVGPAKVIVQLVTNGKNIHLHAHSLVGKHCEDGICTVTAGPKDMVVGFANLGILHVTKKKVFETLEARMTEACIRGYNPGLLVHPDLAYLQAEGGGDRQLGDREKELIRQAALQQTKEMDLSVVRLMFTAFLPDSTGSFTRRLEPVVSDAIYDSKAPNASNLKIVRMDRTAGCVTGGEEIYLLCDKVQKDDIQIRFYEEEENGGVWEGFGDFSPTDVHRQFAIVFKTPKYKDINITKPASVFVQLRRKSDLETSEPKPFLYYPEIKDKEEVQRKRQKLMPNFSDSFGGGSGAGAGGGGMFGSGGGGGGTGSTGPGYSFPHYGFPTYGGITFHPGTTKSNAGMKHGTMDTESKKDPEGCDKSDDKNTVNLFGKVIETTEQDQEPSEATVGNGEVTLTYATGTKEESAGVQDNLFLEKAMQLAKRHANALFDYAVTGDVKMLLAVQRHLTAVQDENGDSVLHLAIIHLHSQLVRDLLEVTSGLISDDIINMRNDLYQTPLHLAVITKQEDVVEDLLRAGADLSLLDRLGNSVLHLAAKEGHDKVLSILLKHKKAALLLDHPNGDGLNAIHLAMMSNSLPCLLLLVAAGADVNAQEQKSGRTALHLAVEHDNISLAGCLLLEGDAHVDSTTYDGTTPLHIAAGRGSTRLAALLKAAGADPLVENFEPLYDLDDSWENAGEDEGVVPGTTPLDMATSWQVFDILNGKPYEPEFTSDDLLAQGDMKQLAEDVKLQLYKLLEIPDPDKNWATLAQKLGLGILNNAFRLSPAPSKTLMDNYEVSGGTVRELVEALRQMGYTEAIEVIQAASSPVKTTSQAHSLPLSPASTRQQIDELRDSDSVCDSGVETSFRKLSFTESLTSGASLLTLNKMPHDYGQEGPLEGKIXXMAEDDPYLGRPEQMFHLDPSLTHTIFNPEVFQPQMALPTDGPYLQILEQPKQRGFRFRYVCEGPSHGGLPGASSEKNKKSYPQVKICNYVGPAKVIVQLVTNGKNIHLHAHSLVGKHCEDGICTVTAGPKDMVVGFANLGILHVTKKKVFETLEARMTEACIRGYNPGLLVHPDLAYLQAEGGGDRQLGDREKELIRQAALQQTKEMDLSVVRLMFTAFLPDSTGSFTRRLEPVVSDAIYDSKAPNASNLKIVRMDRTAGCVTGGEEIYLLCDKVQKDDIQIRFYEEEENGGVWEGFGDFSPTDVHRQFAIVFKTPKYKDINITKPASVFVQLRRKSDLETSEPKPFLYYPEIKDKEEVQRKRQKLMPNFSDSFGGGSGAGAGGGGMFGSGGGGGGTGSTGPGYSFPHYGFPTYGGITFHPGTTKSNAGMKHGTMDTESKKDPEGCDKSDDKNTVNLFGKVIETTEQDQEPSEATVGNGEVTLTYATGTKEESAGVQDNLFLEKAMQLAKRHANALFDYAVTGDVKMLLAVQRHLTAVQDENGDSVLHLAIIHLHSQLVRDLLEVTSGLISDDIINMRNDLYQTPLHLAVITKQEDVVEDLLRAGADLSLLDRLGNSVLHLAAKEGHDKVLSILLKHKKAALLLDHPNGDGLNAIHLAMMSNSLPCLLLLVAAGADVNAQEQKSGRTALHLAVEHDNISLAGCLLLEGDAHVDSTTYDGTTPLHIAAGRGSTRLAALLKAAGADPLVENFEPLYDLDDSWENAGEDEGVVPGTTPLDMATSWQVFDILNGKPYEPEFTSDDLLAQGDMKQLAEDVKLQLYKLLEIPDPDKNWATLAQKLGLGILNNAFRLSPAPSKTLMDNYEVSGGTVRELVEALRQMGYTEAIEVIQAASSPVKTTSQAHSLPLSPASTRQQIDELRDSDSVCDSGVETSFRKLSFTESLTSGASLLTLNKMPHDYGQEGPLEGKIX"
## [6,] "MAEDDPYLGRPEQMFHLDPSLTHTIFNPEVFQPQMALPTADGPYLQILEQPKQRGFRFRYVCEGPSHGGLPGASSEKNKKSYPQVKICNYVGPAKVIVQLVTNGKNIHLHAHSLVGKHCEDGICTVTAGPKDMVVGFANLGILHVTKKKVFETLEARMTEACIRGYNPGLLVHPDLAYLQAEGGGDRQLGDREKELIRQAALQQTKEMDLSVVRLMFTAFLPDSTGSFTRRLEPVVSDAIYDSKAPNASNLKIVRMDRTAGCVTGGEEIYLLCDKVQKDDIQIRFYEEEENGGVWEGFGDFSPTDVHRQFAIVFKTPKYKDINITKPASVFVQLRRKSDLETSEPKPFLYYPEIKDKEEVQRKRQKLMPNFSDSFGGGSGAGAGGGGMFGSGGGGGGTGSTGPGYSFPHYGFPTYGGITFHPGTTKSNAGMKHGTMDTESKKDPEGCDKSDDKNTVNLFGKVIETTEQDQEPSEATVGNGEVTLTYATGTKEESAGVQDNLFLEKAMQLAKRHANALFDYAVTGDVKMLLAVQRHLTAVQDENGDSVLHLAIIHLHSQLVRDLLEVTSGLISDDIINMRNDLYQTPLHLAVITKQEDVVEDLLRAGADLSLLDRLGNSVLHLAAKEGHDKVLSILLKHKKAALLLDHPNGDGLNAIHLAMMSNSLPCLLLLVAAGADVNAQEQKSGRTALHLAVEHDNISLAGCLLLEGDAHVDSTTYDGTTPLHIAAGRGSTRLAALLKAAGADPLVENFEPLYDLDDSWENAGEDEGVVPGTTPLDMATSWQVFDILNGKPYEPEFTSDDLLAQGDMKQLAEDVKLQLYKLLEIPDPDKNWATLAQKLGLGILNNAFRLSPAPSKTLMDNYEVSGGTVRELVEALRQMGYTEAIEVIQAASSPVKTTSQAHSLPLSPASTRQQIDELRDSDSVCDSGVETSFRKLSFTESLTSGASLLTLNKMPHDYGQEGPLEGKIXXMAEDDPYLGRPEQMFHLDPSLTHTIFNPEVFQPQMALPTDGPYLQILEQPKQRGFRFRYVCEGPSHGGLPGASSEKNKKSYPQVKICNYVGPAKVIVQLVTNGKNIHLHAHSLVGKHCEDGICTVTAGPKDMVVGFANLGILHVTKKKVFETLEARMTEACIRGYNPGLLVHPDLAYLQAEGGGDRQLGDREKELIRQAALQQTKEMDLSVVRLMFTAFLPDSTGSFTRRLEPVVSDAIYDSKAPNASNLKIVRMDRTAGCVTGGEEIYLLCDKVQKDDIQIRFYEEEENGGVWEGFGDFSPTDVHRQFAIVFKTPKYKDINITKPASVFVQLRRKSDLETSEPKPFLYYPEIKDKEEVQRKRQKLMPNFSDSFGGGSGAGAGGGGMFGSGGGGGGTGSTGPGYSFPHYGFPTYGGITFHPGTTKSNAGMKHGTMDTESKKDPEGCDKSDDKNTVNLFGKVIETTEQDQEPSEATVGNGEVTLTYATGTKEESAGVQDNLFLEKAMQLAKRHANALFDYAVTGDVKMLLAVQRHLTAVQDENGDSVLHLAIIHLHSQLVRDLLEVTSGLISDDIINMRNDLYQTPLHLAVITKQEDVVEDLLRAGADLSLLDRLGNSVLHLAAKEGHDKVLSILLKHKKAALLLDHPNGDGLNAIHLAMMSNSLPCLLLLVAAGADVNAQEQKSGRTALHLAVEHDNISLAGCLLLEGDAHVDSTTYDGTTPLHIAAGRGSTRLAALLKAAGADPLVENFEPLYDLDDSWENAGEDEGVVPGTTPLDMATSWQVFDILNGKPYEPEFTSDDLLAQGDMKQLAEDVKLQLYKLLEIPDPDKNWATLAQKLGLGILNNAFRLSPAPSKTLMDNYEVSGGTVRELVEALRQMGYTEAIEVIQAASSPVKTTSQAHSLPLSPASTRQQIDELRDSDSVCDSGVETSFRKLSFTESLTSGASLLTLNKMPHDYGQEGPLEGKIX"
##      Mutant_Peptide                Total_RNA Tumor_RNA_Ratio Tumor_RNA
## [1,] "MAEDDPYLGRPEKMFHLDPSLTHTIFN" "NA"      "NA"            "NA"     
## [2,] "MAEDDPYLGRPEKMFHLDPSLTHTIFN" "NA"      "NA"            "NA"     
## [3,] "MAEDDPYLGRPEKMFHLDPSLTHTIFN" "NA"      "NA"            "NA"     
## [4,] "MAEDDPYLGRPEKMFHLDPSLTHTIFN" "NA"      "NA"            "NA"     
## [5,] "MAEDDPYLGRPEKMFHLDPSLTHTIFN" "NA"      "NA"            "NA"     
## [6,] "MAEDDPYLGRPEKMFHLDPSLTHTIFN" "NA"      "NA"            "NA"     
##      Tumor_RNA_based_on_DNA MutRatio MutRatio_Min MutRatio_Max P_I P_R P 
## [1,] "NA"                   "NA"     "NA"         "NA"         NA  NA  NA
## [2,] "NA"                   "NA"     "NA"         "NA"         NA  NA  NA
## [3,] "NA"                   "NA"     "NA"         "NA"         NA  NA  NA
## [4,] "NA"                   "NA"     "NA"         "NA"         NA  NA  NA
## [5,] "NA"                   "NA"     "NA"         "NA"         NA  NA  NA
## [6,] "NA"                   "NA"     "NA"         "NA"         NA  NA  NA
```


``` r
  print(Export_Summary_Fragments(Result_HLA1_Seq, Mut_Rank_th = 5))
```

```
## [[1]]
##                                                     atggcagaag
## Num_Peptide_Per_Grp                                     63.000
## Num_Cond_Peptide_Per_Grp                                63.000
## Num_Rest_Peptide_Per_Grp                                15.000
## Num_Rest_Peptide_Per_Grp / Num_Cond_Peptide_Per_Grp      0.238
## -logP                                                    0.621
## 
## [[2]]
##                                                         
## Num_Peptide_Per_NM                                63.000
## Num_Cond_Peptide_Per_NM                           63.000
## Num_Rest_Peptide_Per_NM                           15.000
## Num_Rest_Peptide_Per_NM / Num_Cond_Peptide_Per_NM  0.238
## -logP                                              0.621
## 
## [[3]]
##                                                     MAEDDPYLGRPEKMFHLDPSLTHTIFN-0_atggcagaag-0_1
## Num_Peptide_Per_Pep                                                                       63.000
## Num_Cond_Peptide_Per_Pep                                                                  63.000
## Num_Rest_Peptide_Per_Pep                                                                  15.000
## Num_Rest_Peptide_Per_Pep / Num_Cond_Peptide_Per_Pep                                        0.238
## -logP                                                                                      0.621
```

To calculate the binding affinity of neoantigen candaites, which are directly generated from RNA sequences, to HLA ClassII. The peptides included in the riginal genes ("NFKB1", "BCL3") are removed from the results.


``` r
  Result_HLA2_Seq <- MainSeqFragmentClass2(input_sequence = "atggcagaagatgatccatatttgggaaggcctgaaaaaatgtttcatttggatccttctttgactcatacaatatttaatc",
                                           file_name_in_hla_table = "sample",
                                           hla_file = "./../data/sample_hla_table_c2.txt",
                                           hmdir = getwd(),
                                           job_id = "NO_job_id",
                                           refflat_file  = "./../lib/refFlat.grch37.txt",
                                           refmrna_file = "./../lib/refMrna.grch37.fa",
                                           netMHCIIpan_dir = "./../lib/netMHCIIpan-4.3/netMHCIIpan",
                                           reference_gene_symbol = c("NFKB1", "BCL3"))
```

```
## [1] "Sequence: atggcagaagatgatccatatttgggaaggcctgaaaaaatgtttcatttggatccttctttgactcatacaatatttaatc"
##  [1] "HLAtype: DPA1*01:03" "HLAtype: DPA1*02:01" "HLAtype: DPB1*02:01"
##  [4] "HLAtype: DPB1*09:01" "HLAtype: DQA1*01:02" "HLAtype: DQA1*05:05"
##  [7] "HLAtype: DQB1*03:01" "HLAtype: DQB1*06:04" "HLAtype: DRB1*11:04"
## [10] "HLAtype: DRB1*13:02"
## [1] "refFLAT: ./../lib/refFlat.grch37.txt"
## [1] "refMrna: ./../lib/refMrna.grch37.fa"
## [1] "Wt-NM_ID:  NA"
## [1] "Wt-Gene Symbol: NFKB1" "Wt-Gene Symbol: BCL3" 
## [1] "Start Position of Reading Frame is: 1"
## [1] "Peptide Successfully Generated!!"
## [1] "Executing netMHCpan to result.NO_job_id.SeqFragment2"
## [1] "./../lib/netMHCIIpan-4.3/netMHCIIpan -length 15 -f result.NO_job_id.SeqFragment2/NO_job_id_SeqFragment.peptide.fasta -choose -cha DPA10103 -choose -chb DPB10201 > result.NO_job_id.SeqFragment2/NA.NO_job_id_SeqFragment.HLACLASS2.1.peptide.txt"
## [1] "./../lib/netMHCIIpan-4.3/netMHCIIpan -length 15 -f result.NO_job_id.SeqFragment2/NO_job_id_SeqFragment.peptide.fasta -choose -cha DPA10103 -choose -chb DPB10901 > result.NO_job_id.SeqFragment2/NA.NO_job_id_SeqFragment.HLACLASS2.2.peptide.txt"
## [1] "./../lib/netMHCIIpan-4.3/netMHCIIpan -length 15 -f result.NO_job_id.SeqFragment2/NO_job_id_SeqFragment.peptide.fasta -choose -cha DPA10201 -choose -chb DPB10201 > result.NO_job_id.SeqFragment2/NA.NO_job_id_SeqFragment.HLACLASS2.3.peptide.txt"
## [1] "./../lib/netMHCIIpan-4.3/netMHCIIpan -length 15 -f result.NO_job_id.SeqFragment2/NO_job_id_SeqFragment.peptide.fasta -choose -cha DPA10201 -choose -chb DPB10901 > result.NO_job_id.SeqFragment2/NA.NO_job_id_SeqFragment.HLACLASS2.4.peptide.txt"
## [1] "./../lib/netMHCIIpan-4.3/netMHCIIpan -length 15 -f result.NO_job_id.SeqFragment2/NO_job_id_SeqFragment.peptide.fasta -choose -cha DQA10102 -choose -chb DQB10301 > result.NO_job_id.SeqFragment2/NO_job_id_SeqFragment.HLACLASS2.5.peptide.txt"
## [1] "./../lib/netMHCIIpan-4.3/netMHCIIpan -length 15 -f result.NO_job_id.SeqFragment2/NO_job_id_SeqFragment.peptide.fasta -choose -cha DQA10102 -choose -chb DQB10604 > result.NO_job_id.SeqFragment2/NO_job_id_SeqFragment.HLACLASS2.6.peptide.txt"
## [1] "./../lib/netMHCIIpan-4.3/netMHCIIpan -length 15 -f result.NO_job_id.SeqFragment2/NO_job_id_SeqFragment.peptide.fasta -choose -cha DQA10505 -choose -chb DQB10301 > result.NO_job_id.SeqFragment2/NO_job_id_SeqFragment.HLACLASS2.7.peptide.txt"
## [1] "./../lib/netMHCIIpan-4.3/netMHCIIpan -length 15 -f result.NO_job_id.SeqFragment2/NO_job_id_SeqFragment.peptide.fasta -choose -cha DQA10505 -choose -chb DQB10604 > result.NO_job_id.SeqFragment2/NO_job_id_SeqFragment.HLACLASS2.8.peptide.txt"
## [1] "./../lib/netMHCIIpan-4.3/netMHCIIpan -length 15 -f result.NO_job_id.SeqFragment2/NO_job_id_SeqFragment.peptide.fasta -a DRB1_1104 > result.NO_job_id.SeqFragment2/NA.NO_job_id_SeqFragment.HLACLASS2.9.peptide.txt"
## [1] "./../lib/netMHCIIpan-4.3/netMHCIIpan -length 15 -f result.NO_job_id.SeqFragment2/NO_job_id_SeqFragment.peptide.fasta -a DRB1_1302 > result.NO_job_id.SeqFragment2/NA.NO_job_id_SeqFragment.HLACLASS2.10.peptide.txt"
## [1] "Merging Results..."
## [1] "/Users/takanorihasegawa/Documents/GitHub/Neoantimon/vignettes/result.NO_job_id.SeqFragment2/NA.NO_job_id_SeqFragment.HLACLASS2.1.peptide.txt"
## [1] "100 perc. fin"
## [1] "/Users/takanorihasegawa/Documents/GitHub/Neoantimon/vignettes/result.NO_job_id.SeqFragment2/NA.NO_job_id_SeqFragment.HLACLASS2.10.peptide.txt"
## [1] "100 perc. fin"
## [1] "/Users/takanorihasegawa/Documents/GitHub/Neoantimon/vignettes/result.NO_job_id.SeqFragment2/NA.NO_job_id_SeqFragment.HLACLASS2.2.peptide.txt"
## [1] "100 perc. fin"
## [1] "/Users/takanorihasegawa/Documents/GitHub/Neoantimon/vignettes/result.NO_job_id.SeqFragment2/NA.NO_job_id_SeqFragment.HLACLASS2.3.peptide.txt"
## [1] "100 perc. fin"
## [1] "/Users/takanorihasegawa/Documents/GitHub/Neoantimon/vignettes/result.NO_job_id.SeqFragment2/NA.NO_job_id_SeqFragment.HLACLASS2.4.peptide.txt"
## [1] "100 perc. fin"
## [1] "/Users/takanorihasegawa/Documents/GitHub/Neoantimon/vignettes/result.NO_job_id.SeqFragment2/NA.NO_job_id_SeqFragment.HLACLASS2.9.peptide.txt"
## [1] "100 perc. fin"
## [1] "/Users/takanorihasegawa/Documents/GitHub/Neoantimon/vignettes/result.NO_job_id.SeqFragment2/NO_job_id_SeqFragment.HLACLASS2.5.peptide.txt"
## [1] "100 perc. fin"
## [1] "/Users/takanorihasegawa/Documents/GitHub/Neoantimon/vignettes/result.NO_job_id.SeqFragment2/NO_job_id_SeqFragment.HLACLASS2.6.peptide.txt"
## [1] "100 perc. fin"
## [1] "/Users/takanorihasegawa/Documents/GitHub/Neoantimon/vignettes/result.NO_job_id.SeqFragment2/NO_job_id_SeqFragment.HLACLASS2.7.peptide.txt"
## [1] "100 perc. fin"
## [1] "/Users/takanorihasegawa/Documents/GitHub/Neoantimon/vignettes/result.NO_job_id.SeqFragment2/NO_job_id_SeqFragment.HLACLASS2.8.peptide.txt"
## [1] "100 perc. fin"
## [1] "Successfully Finished."
```

``` r
  Result_HLA2_Seq <- CalculatePriorityScores(result = Result_HLA2_Seq, useRNAvaf = FALSE)
```

```
## Warning in CalculatePriorityScores(result = Result_HLA2_Seq, useRNAvaf =
## FALSE): NAs introduced by coercion
```

``` r
  print(head(Result_HLA2_Seq))
```

```
##      HLA         Pos Gene           Evaluated_Mutant_Peptide_Core
## [1,] "DRB1_0101" "1" "0_atggcagaag" "PYLGRPEKM"                  
## [2,] "DRB1_0101" "2" "0_atggcagaag" "YLGRPEKMF"                  
## [3,] "DRB1_0101" "3" "0_atggcagaag" "YLGRPEKMF"                  
## [4,] "DRB1_0101" "4" "0_atggcagaag" "YLGRPEKMF"                  
## [5,] "DRB1_0101" "5" "0_atggcagaag" "YLGRPEKMF"                  
## [6,] "DRB1_0101" "6" "0_atggcagaag" "GRPEKMFHL"                  
##      Evaluated_Mutant_Peptide Mut_EL     Mut_Rank Chr NM_ID ReadingFrame
## [1,] "MAEDDPYLGRPEKMF"        "0.000116" "74.92"  "0" ""    "1"         
## [2,] "AEDDPYLGRPEKMFH"        "0.000639" "53.83"  "0" ""    "1"         
## [3,] "EDDPYLGRPEKMFHL"        "0.000790" "51.22"  "0" ""    "1"         
## [4,] "DDPYLGRPEKMFHLD"        "0.001267" "45.36"  "0" ""    "1"         
## [5,] "DPYLGRPEKMFHLDP"        "0.000216" "67.70"  "0" ""    "1"         
## [6,] "PYLGRPEKMFHLDPS"        "0.000048" "84.48"  "0" ""    "1"         
##      SequenceNumber Chrs                                      
## [1,] "1"            "chr19;chr4;chr4;chr4;chr4;chr4;chr4;chr4"
## [2,] "1"            "chr19;chr4;chr4;chr4;chr4;chr4;chr4;chr4"
## [3,] "1"            "chr19;chr4;chr4;chr4;chr4;chr4;chr4;chr4"
## [4,] "1"            "chr19;chr4;chr4;chr4;chr4;chr4;chr4;chr4"
## [5,] "1"            "chr19;chr4;chr4;chr4;chr4;chr4;chr4;chr4"
## [6,] "1"            "chr19;chr4;chr4;chr4;chr4;chr4;chr4;chr4"
##      NM_IDs                                                                                             
## [1,] "NM_005178;NM_001319226;NM_003998;NM_001165412;NM_001382626;NM_001382627;NM_001382628;NM_001382625"
## [2,] "NM_005178;NM_001319226;NM_003998;NM_001165412;NM_001382626;NM_001382627;NM_001382628;NM_001382625"
## [3,] "NM_005178;NM_001319226;NM_003998;NM_001165412;NM_001382626;NM_001382627;NM_001382628;NM_001382625"
## [4,] "NM_005178;NM_001319226;NM_003998;NM_001165412;NM_001382626;NM_001382627;NM_001382628;NM_001382625"
## [5,] "NM_005178;NM_001319226;NM_003998;NM_001165412;NM_001382626;NM_001382627;NM_001382628;NM_001382625"
## [6,] "NM_005178;NM_001319226;NM_003998;NM_001165412;NM_001382626;NM_001382627;NM_001382628;NM_001382625"
##      GeneIDs                                         
## [1,] "BCL3;NFKB1;NFKB1;NFKB1;NFKB1;NFKB1;NFKB1;NFKB1"
## [2,] "BCL3;NFKB1;NFKB1;NFKB1;NFKB1;NFKB1;NFKB1;NFKB1"
## [3,] "BCL3;NFKB1;NFKB1;NFKB1;NFKB1;NFKB1;NFKB1;NFKB1"
## [4,] "BCL3;NFKB1;NFKB1;NFKB1;NFKB1;NFKB1;NFKB1;NFKB1"
## [5,] "BCL3;NFKB1;NFKB1;NFKB1;NFKB1;NFKB1;NFKB1;NFKB1"
## [6,] "BCL3;NFKB1;NFKB1;NFKB1;NFKB1;NFKB1;NFKB1;NFKB1"
##      Exon_Starts                                                                     
## [1,] "45251964;103423090;103422515;103422515;103422515;103422515;103422515;103423090"
## [2,] "45251964;103423090;103422515;103422515;103422515;103422515;103422515;103423090"
## [3,] "45251964;103423090;103422515;103422515;103422515;103422515;103422515;103423090"
## [4,] "45251964;103423090;103422515;103422515;103422515;103422515;103422515;103423090"
## [5,] "45251964;103423090;103422515;103422515;103422515;103422515;103422515;103423090"
## [6,] "45251964;103423090;103422515;103422515;103422515;103422515;103422515;103423090"
##      Exon_Ends                                                                       
## [1,] "45263301;103538459;103538459;103538459;103538459;103538459;103538459;103538459"
## [2,] "45263301;103538459;103538459;103538459;103538459;103538459;103538459;103538459"
## [3,] "45263301;103538459;103538459;103538459;103538459;103538459;103538459;103538459"
## [4,] "45263301;103538459;103538459;103538459;103538459;103538459;103538459;103538459"
## [5,] "45263301;103538459;103538459;103538459;103538459;103538459;103538459;103538459"
## [6,] "45263301;103538459;103538459;103538459;103538459;103538459;103538459;103538459"
##      GroupID NumOfPeptides NumOfStops
## [1,] "0_1"   "27"          "0"       
## [2,] "0_1"   "27"          "0"       
## [3,] "0_1"   "27"          "0"       
## [4,] "0_1"   "27"          "0"       
## [5,] "0_1"   "27"          "0"       
## [6,] "0_1"   "27"          "0"       
##      Wt_Peptide                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           
## [1,] "MPRCPAGAMDEGPVDLRTRPKAAGLPGAALPLRKRPLRAPSPEPAAPRGAAGLVVPLDPLRGGCDLPAVPGPPHGLARPEALYYPGALLPLYPTRAMGSPFPLVNLPTPLYPMMCPMEHPLSADIAMATRADEDGDTPLHIAVVQGNLPAVHRLVNLFQQGGRELDIYNNLRQTPLHLAVITTLPSVVRLLVTAGASPMALDRHGQTAAHLACEHRSPTCLRALLDSAAPGTLDLEARNYDGLTALHVAVNTECQETVQLLLERGADIDAVDIKSGRSPLIHAVENNSLSMVQLLLQHGANVNAQMYSGSSALHSASGRGLLPLVRTLVRSGADSSLKNCHNDTPLMVARSRRVIDILRGKATRPASTSQPDPSPDRSANTSPESSSRLSSNGLLSASPSSSPSQSPPRDPPGFPMAPPNFFLPSPSPPAFLPFAGVLRGPGRPVPPSPAPGGSXXMAEDDPYLGRPEQMFHLDPSLTHTIFNPEVFQPQMALPTDGPYLQILEQPKQRGFRFRYVCEGPSHGGLPGASSEKNKKSYPQVKICNYVGPAKVIVQLVTNGKNIHLHAHSLVGKHCEDGICTVTAGPKDMVVGFANLGILHVTKKKVFETLEARMTEACIRGYNPGLLVHPDLAYLQAEGGGDRQLGDREKELIRQAALQQTKEMDLSVVRLMFTAFLPDSTGSFTRRLEPVVSDAIYDSKAPNASNLKIVRMDRTAGCVTGGEEIYLLCDKVQKDDIQIRFYEEEENGGVWEGFGDFSPTDVHRQFAIVFKTPKYKDINITKPASVFVQLRRKSDLETSEPKPFLYYPEIKDKEEVQRKRQKLMPNFSDSFGGGSGAGAGGGGMFGSGGGGGGTGSTGPGYSFPHYGFPTYGGITFHPGTTKSNAGMKHGTMDTESKKDPEGCDKSDDKNTVNLFGKVIETTEQDQEPSEATVGNGEVTLTYATGTKEESAGVQDNLFLEKAMQLAKRHANALFDYAVTGDVKMLLAVQRHLTAVQDENGDSVLHLAIIHLHSQLVRDLLEVTSGLISDDIINMRNDLYQTPLHLAVITKQEDVVEDLLRAGADLSLLDRLGNSVLHLAAKEGHDKVLSILLKHKKAALLLDHPNGDGLNAIHLAMMSNSLPCLLLLVAAGADVNAQEQKSGRTALHLAVEHDNISLAGCLLLEGDAHVDSTTYDGTTPLHIAAGRGSTRLAALLKAAGADPLVENFEPLYDLDDSWENAGEDEGVVPGTTPLDMATSWQVFDILNGKPYEPEFTSDDLLAQGDMKQLAEDVKLQLYKLLEIPDPDKNWATLAQKLGLGILNNAFRLSPAPSKTLMDNYEVSGGTVRELVEALRQMGYTEAIEVIQAASSPVKTTSQAHSLPLSPASTRQQIDELRDSDSVCDSGVETSFRKLSFTESLTSGASLLTLNKMPHDYGQEGPLEGKIXXMAEDDPYLGRPEQMFHLDPSLTHTIFNPEVFQPQMALPTADGPYLQILEQPKQRGFRFRYVCEGPSHGGLPGASSEKNKKSYPQVKICNYVGPAKVIVQLVTNGKNIHLHAHSLVGKHCEDGICTVTAGPKDMVVGFANLGILHVTKKKVFETLEARMTEACIRGYNPGLLVHPDLAYLQAEGGGDRQLGDREKELIRQAALQQTKEMDLSVVRLMFTAFLPDSTGSFTRRLEPVVSDAIYDSKAPNASNLKIVRMDRTAGCVTGGEEIYLLCDKVQKDDIQIRFYEEEENGGVWEGFGDFSPTDVHRQFAIVFKTPKYKDINITKPASVFVQLRRKSDLETSEPKPFLYYPEIKDKEEVQRKRQKLMPNFSDSFGGGSGAGAGGGGMFGSGGGGGGTGSTGPGYSFPHYGFPTYGGITFHPGTTKSNAGMKHGTMDTESKKDPEGCDKSDDKNTVNLFGKVIETTEQDQEPSEATVGNGEVTLTYATGTKEESAGVQDNLFLEKAMQLAKRHANALFDYAVTGDVKMLLAVQRHLTAVQDENGDSVLHLAIIHLHSQLVRDLLEVTSGLISDDIINMRNDLYQTPLHLAVITKQEDVVEDLLRAGADLSLLDRLGNSVLHLAAKEGHDKVLSILLKHKKAALLLDHPNGDGLNAIHLAMMSNSLPCLLLLVAAGADVNAQEQKSGRTALHLAVEHDNISLAGCLLLEGDAHVDSTTYDGTTPLHIAAGRGSTRLAALLKAAGADPLVENFEPLYDLDDSWENAGEDEGVVPGTTPLDMATSWQVFDILNGKPYEPEFTSDDLLAQGDMKQLAEDVKLQLYKLLEIPDPDKNWATLAQKLGLGILNNAFRLSPAPSKTLMDNYEVSGGTVRELVEALRQMGYTEAIEVIQAASSPVKTTSQAHSLPLSPASTRQQIDELRDSDSVCDSGVETSFRKLSFTESLTSGASLLTLNKMPHDYGQEGPLEGKIXXMAEDDPYLGRPEQMFHLDPSLTHTIFNPEVFQPQMALPTDGPYLQILEQPKQRGFRFRYVCEGPSHGGLPGASSEKNKKSYPQVKICNYVGPAKVIVQLVTNGKNIHLHAHSLVGKHCEDGICTVTAGPKDMVVGFANLGILHVTKKKVFETLEARMTEACIRGYNPGLLVHPDLAYLQAEGGGDRQLGDREKELIRQAALQQTKEMDLSVVRLMFTAFLPDSTGSFTRRLEPVVSDAIYDSKAPNASNLKIVRMDRTAGCVTGGEEIYLLCDKVQKDDIQIRFYEEEENGGVWEGFGDFSPTDVHRQFAIVFKTPKYKDINITKPASVFVQLRRKSDLETSEPKPFLYYPEIKDKEEVQRKRQKLMPNFSDSFGGGSGAGAGGGGMFGSGGGGGGTGSTGPGYSFPHYGFPTYGGITFHPGTTKSNAGMKHGTMDTESKKDPEGCDKSDDKNTVNLFGKVIETTEQDQEPSEATVGNGEVTLTYATGTKEESAGVQDNLFLEKAMQLAKRHANALFDYAVTGDVKMLLAVQRHLTAVQDENGDSVLHLAIIHLHSQLVRDLLEVTSGLISDDIINMRNDLYQTPLHLAVITKQEDVVEDLLRAGADLSLLDRLGNSVLHLAAKEGHDKVLSILLKHKKAALLLDHPNGDGLNAIHLAMMSNSLPCLLLLVAAGADVNAQEQKSGRTALHLAVEHDNISLAGCLLLEGDAHVDSTTYDGTTPLHIAAGRGSTRLAALLKAAGADPLVENFEPLYDLDDSWENAGEDEGVVPGTTPLDMATSWQVFDILNGKPYEPEFTSDDLLAQGDMKQLAEDVKLQLYKLLEIPDPDKNWATLAQKLGLGILNNAFRLSPAPSKTLMDNYEVSGGTVRELVEALRQMGYTEAIEVIQAASSPVKTTSQAHSLPLSPASTRQQIDELRDSDSVCDSGVETSFRKLSFTESLTSGASLLTLNKMPHDYGQEGPLEGKIXXMAEDDPYLGRPEQMFHLDPSLTHTIFNPEVFQPQMALPTADGPYLQILEQPKQRGFRFRYVCEGPSHGGLPGASSEKNKKSYPQVKICNYVGPAKVIVQLVTNGKNIHLHAHSLVGKHCEDGICTVTAGPKDMVVGFANLGILHVTKKKVFETLEARMTEACIRGYNPGLLVHPDLAYLQAEGGGDRQLGDREKELIRQAALQQTKEMDLSVVRLMFTAFLPDSTGSFTRRLEPVVSDAIYDSKAPNASNLKIVRMDRTAGCVTGGEEIYLLCDKVQKDDIQIRFYEEEENGGVWEGFGDFSPTDVHRQFAIVFKTPKYKDINITKPASVFVQLRRKSDLETSEPKPFLYYPEIKDKEEVQRKRQKLMPNFSDSFGGGSGAGAGGGGMFGSGGGGGGTGSTGPGYSFPHYGFPTYGGITFHPGTTKSNAGMKHGTMDTESKKDPEGCDKSDDKNTVNLFGKVIETTEQDQEPSEATVGNGEVTLTYATGTKEESAGVQDNLFLEKAMQLAKRHANALFDYAVTGDVKMLLAVQRHLTAVQDENGDSVLHLAIIHLHSQLVRDLLEVTSGLISDDIINMRNDLYQTPLHLAVITKQEDVVEDLLRAGADLSLLDRLGNSVLHLAAKEGHDKVLSILLKHKKAALLLDHPNGDGLNAIHLAMMSNSLPCLLLLVAAGADVNAQEQKSGRTALHLAVEHDNISLAGCLLLEGDAHVDSTTYDGTTPLHIAAGRGSTRLAALLKAAGADPLVENFEPLYDLDDSWENAGEDEGVVPGTTPLDMATSWQVFDILNGKPYEPEFTSDDLLAQGDMKQLAEDVKLQLYKLLEIPDPDKNWATLAQKLGLGILNNAFRLSPAPSKTLMDNYEVSGGTVRELVEALRQMGYTEAIEVIQAASSPVKTTSQAHSLPLSPASTRQQIDELRDSDSVCDSGVETSFRKLSFTESLTSGASLLTLNKMPHDYGQEGPLEGKIXXMAEDDPYLGRPEQMFHLDPSLTHTIFNPEVFQPQMALPTDGPYLQILEQPKQRGFRFRYVCEGPSHGGLPGASSEKNKKSYPQVKICNYVGPAKVIVQLVTNGKNIHLHAHSLVGKHCEDGICTVTAGPKDMVVGFANLGILHVTKKKVFETLEARMTEACIRGYNPGLLVHPDLAYLQAEGGGDRQLGDREKELIRQAALQQTKEMDLSVVRLMFTAFLPDSTGSFTRRLEPVVSDAIYDSKAPNASNLKIVRMDRTAGCVTGGEEIYLLCDKVQKDDIQIRFYEEEENGGVWEGFGDFSPTDVHRQFAIVFKTPKYKDINITKPASVFVQLRRKSDLETSEPKPFLYYPEIKDKEEVQRKRQKLMPNFSDSFGGGSGAGAGGGGMFGSGGGGGGTGSTGPGYSFPHYGFPTYGGITFHPGTTKSNAGMKHGTMDTESKKDPEGCDKSDDKNTVNLFGKVIETTEQDQEPSEATVGNGEVTLTYATGTKEESAGVQDNLFLEKAMQLAKRHANALFDYAVTGDVKMLLAVQRHLTAVQDENGDSVLHLAIIHLHSQLVRDLLEVTSGLISDDIINMRNDLYQTPLHLAVITKQEDVVEDLLRAGADLSLLDRLGNSVLHLAAKEGHDKVLSILLKHKKAALLLDHPNGDGLNAIHLAMMSNSLPCLLLLVAAGADVNAQEQKSGRTALHLAVEHDNISLAGCLLLEGDAHVDSTTYDGTTPLHIAAGRGSTRLAALLKAAGADPLVENFEPLYDLDDSWENAGEDEGVVPGTTPLDMATSWQVFDILNGKPYEPEFTSDDLLAQGDMKQLAEDVKLQLYKLLEIPDPDKNWATLAQKLGLGILNNAFRLSPAPSKTLMDNYEVSGGTVRELVEALRQMGYTEAIEVIQAASSPVKTTSQAHSLPLSPASTRQQIDELRDSDSVCDSGVETSFRKLSFTESLTSGASLLTLNKMPHDYGQEGPLEGKIXXMFHLDPSLTHTIFNPEVFQPQMALPTDGPYLQILEQPKQRGFRFRYVCEGPSHGGLPGASSEKNKKSYPQVKICNYVGPAKVIVQLVTNGKNIHLHAHSLVGKHCEDGICTVTAGPKDMVVGFANLGILHVTKKKVFETLEARMTEACIRGYNPGLLVHPDLAYLQAEGGGDRQLGDREKELIRQAALQQTKEMDLSVVRLMFTAFLPDSTGSFTRRLEPVVSDAIYDSKAPNASNLKIVRMDRTAGCVTGGEEIYLLCDKVQKDDIQIRFYEEEENGGVWEGFGDFSPTDVHRQFAIVFKTPKYKDINITKPASVFVQLRRKSDLETSEPKPFLYYPEIKDKEEVQRKRQKLMPNFSDSFGGGSGAGAGGGGMFGSGGGGGGTGSTGPGYSFPHYGFPTYGGITFHPGTTKSNAGMKHGTMDTESKKDPEGCDKSDDKNTVNLFGKVIETTEQDQEPSEATVGNGEVTLTYATGTKEESAGVQDNLFLEKAMQLAKRHANALFDYAVTGDVKMLLAVQRHLTAVQDENGDSVLHLAIIHLHSQLVRDLLEVTSGLISDDIINMRNDLYQTPLHLAVITKQEDVVEDLLRAGADLSLLDRLGNSVLHLAAKEGHDKVLSILLKHKKAALLLDHPNGDGLNAIHLAMMSNSLPCLLLLVAAGADVNAQEQKSGRTALHLAVEHDNISLAGCLLLEGDAHVDSTTYDGTTPLHIAAGRGSTRLAALLKAAGADPLVENFEPLYDLDDSWENAGEDEGVVPGTTPLDMATSWQVFDILNGKPYEPEFTSDDLLAQGDMKQLAEDVKLQLYKLLEIPDPDKNWATLAQKLGLGILNNAFRLSPAPSKTLMDNYEVSGGTVRELVEALRQMGYTEAIEVIQAASSPVKTTSQAHSLPLSPASTRQQIDELRDSDSVCDSGVETSFRKLSFTESLTSGASLLTLNKMPHDYGQEGPLEGKIXXMAEDDPYLGRPEQMFHLDPSLTHTIFNPEVFQPQMALPTADGPYLQILEQPKQRGFRFRYVCEGPSHGGLPGASSEKNKKSYPQVKICNYVGPAKVIVQLVTNGKNIHLHAHSLVGKHCEDGICTVTAGPKDMVVGFANLGILHVTKKKVFETLEARMTEACIRGYNPGLLVHPDLAYLQAEGGGDRQLGDREKELIRQAALQQTKEMDLSVVRLMFTAFLPDSTGSFTRRLEPVVSDAIYDSKAPNASNLKIVRMDRTAGCVTGGEEIYLLCDKVQKDDIQIRFYEEEENGGVWEGFGDFSPTDVHRQFAIVFKTPKYKDINITKPASVFVQLRRKSDLETSEPKPFLYYPEIKDKEEVQRKRQKLMPNFSDSFGGGSGAGAGGGGMFGSGGGGGGTGSTGPGYSFPHYGFPTYGGITFHPGTTKSNAGMKHGTMDTESKKDPEGCDKSDDKNTVNLFGKVIETTEQDQEPSEATVGNGEVTLTYATGTKEESAGVQDNLFLEKAMQLAKRHANALFDYAVTGDVKMLLAVQRHLTAVQDENGDSVLHLAIIHLHSQLVRDLLEVTSGLISDDIINMRNDLYQTPLHLAVITKQEDVVEDLLRAGADLSLLDRLGNSVLHLAAKEGHDKVLSILLKHKKAALLLDHPNGDGLNAIHLAMMSNSLPCLLLLVAAGADVNAQEQKSGRTALHLAVEHDNISLAGCLLLEGDAHVDSTTYDGTTPLHIAAGRGSTRLAALLKAAGADPLVENFEPLYDLDDSWENAGEDEGVVPGTTPLDMATSWQVFDILNGKPYEPEFTSDDLLAQGDMKQLAEDVKLQLYKLLEIPDPDKNWATLAQKLGLGILNNAFRLSPAPSKTLMDNYEVSGGTVRELVEALRQMGYTEAIEVIQAASSPVKTTSQAHSLPLSPASTRQQIDELRDSDSVCDSGVETSFRKLSFTESLTSGASLLTLNKMPHDYGQEGPLEGKIX"
## [2,] "MPRCPAGAMDEGPVDLRTRPKAAGLPGAALPLRKRPLRAPSPEPAAPRGAAGLVVPLDPLRGGCDLPAVPGPPHGLARPEALYYPGALLPLYPTRAMGSPFPLVNLPTPLYPMMCPMEHPLSADIAMATRADEDGDTPLHIAVVQGNLPAVHRLVNLFQQGGRELDIYNNLRQTPLHLAVITTLPSVVRLLVTAGASPMALDRHGQTAAHLACEHRSPTCLRALLDSAAPGTLDLEARNYDGLTALHVAVNTECQETVQLLLERGADIDAVDIKSGRSPLIHAVENNSLSMVQLLLQHGANVNAQMYSGSSALHSASGRGLLPLVRTLVRSGADSSLKNCHNDTPLMVARSRRVIDILRGKATRPASTSQPDPSPDRSANTSPESSSRLSSNGLLSASPSSSPSQSPPRDPPGFPMAPPNFFLPSPSPPAFLPFAGVLRGPGRPVPPSPAPGGSXXMAEDDPYLGRPEQMFHLDPSLTHTIFNPEVFQPQMALPTDGPYLQILEQPKQRGFRFRYVCEGPSHGGLPGASSEKNKKSYPQVKICNYVGPAKVIVQLVTNGKNIHLHAHSLVGKHCEDGICTVTAGPKDMVVGFANLGILHVTKKKVFETLEARMTEACIRGYNPGLLVHPDLAYLQAEGGGDRQLGDREKELIRQAALQQTKEMDLSVVRLMFTAFLPDSTGSFTRRLEPVVSDAIYDSKAPNASNLKIVRMDRTAGCVTGGEEIYLLCDKVQKDDIQIRFYEEEENGGVWEGFGDFSPTDVHRQFAIVFKTPKYKDINITKPASVFVQLRRKSDLETSEPKPFLYYPEIKDKEEVQRKRQKLMPNFSDSFGGGSGAGAGGGGMFGSGGGGGGTGSTGPGYSFPHYGFPTYGGITFHPGTTKSNAGMKHGTMDTESKKDPEGCDKSDDKNTVNLFGKVIETTEQDQEPSEATVGNGEVTLTYATGTKEESAGVQDNLFLEKAMQLAKRHANALFDYAVTGDVKMLLAVQRHLTAVQDENGDSVLHLAIIHLHSQLVRDLLEVTSGLISDDIINMRNDLYQTPLHLAVITKQEDVVEDLLRAGADLSLLDRLGNSVLHLAAKEGHDKVLSILLKHKKAALLLDHPNGDGLNAIHLAMMSNSLPCLLLLVAAGADVNAQEQKSGRTALHLAVEHDNISLAGCLLLEGDAHVDSTTYDGTTPLHIAAGRGSTRLAALLKAAGADPLVENFEPLYDLDDSWENAGEDEGVVPGTTPLDMATSWQVFDILNGKPYEPEFTSDDLLAQGDMKQLAEDVKLQLYKLLEIPDPDKNWATLAQKLGLGILNNAFRLSPAPSKTLMDNYEVSGGTVRELVEALRQMGYTEAIEVIQAASSPVKTTSQAHSLPLSPASTRQQIDELRDSDSVCDSGVETSFRKLSFTESLTSGASLLTLNKMPHDYGQEGPLEGKIXXMAEDDPYLGRPEQMFHLDPSLTHTIFNPEVFQPQMALPTADGPYLQILEQPKQRGFRFRYVCEGPSHGGLPGASSEKNKKSYPQVKICNYVGPAKVIVQLVTNGKNIHLHAHSLVGKHCEDGICTVTAGPKDMVVGFANLGILHVTKKKVFETLEARMTEACIRGYNPGLLVHPDLAYLQAEGGGDRQLGDREKELIRQAALQQTKEMDLSVVRLMFTAFLPDSTGSFTRRLEPVVSDAIYDSKAPNASNLKIVRMDRTAGCVTGGEEIYLLCDKVQKDDIQIRFYEEEENGGVWEGFGDFSPTDVHRQFAIVFKTPKYKDINITKPASVFVQLRRKSDLETSEPKPFLYYPEIKDKEEVQRKRQKLMPNFSDSFGGGSGAGAGGGGMFGSGGGGGGTGSTGPGYSFPHYGFPTYGGITFHPGTTKSNAGMKHGTMDTESKKDPEGCDKSDDKNTVNLFGKVIETTEQDQEPSEATVGNGEVTLTYATGTKEESAGVQDNLFLEKAMQLAKRHANALFDYAVTGDVKMLLAVQRHLTAVQDENGDSVLHLAIIHLHSQLVRDLLEVTSGLISDDIINMRNDLYQTPLHLAVITKQEDVVEDLLRAGADLSLLDRLGNSVLHLAAKEGHDKVLSILLKHKKAALLLDHPNGDGLNAIHLAMMSNSLPCLLLLVAAGADVNAQEQKSGRTALHLAVEHDNISLAGCLLLEGDAHVDSTTYDGTTPLHIAAGRGSTRLAALLKAAGADPLVENFEPLYDLDDSWENAGEDEGVVPGTTPLDMATSWQVFDILNGKPYEPEFTSDDLLAQGDMKQLAEDVKLQLYKLLEIPDPDKNWATLAQKLGLGILNNAFRLSPAPSKTLMDNYEVSGGTVRELVEALRQMGYTEAIEVIQAASSPVKTTSQAHSLPLSPASTRQQIDELRDSDSVCDSGVETSFRKLSFTESLTSGASLLTLNKMPHDYGQEGPLEGKIXXMAEDDPYLGRPEQMFHLDPSLTHTIFNPEVFQPQMALPTDGPYLQILEQPKQRGFRFRYVCEGPSHGGLPGASSEKNKKSYPQVKICNYVGPAKVIVQLVTNGKNIHLHAHSLVGKHCEDGICTVTAGPKDMVVGFANLGILHVTKKKVFETLEARMTEACIRGYNPGLLVHPDLAYLQAEGGGDRQLGDREKELIRQAALQQTKEMDLSVVRLMFTAFLPDSTGSFTRRLEPVVSDAIYDSKAPNASNLKIVRMDRTAGCVTGGEEIYLLCDKVQKDDIQIRFYEEEENGGVWEGFGDFSPTDVHRQFAIVFKTPKYKDINITKPASVFVQLRRKSDLETSEPKPFLYYPEIKDKEEVQRKRQKLMPNFSDSFGGGSGAGAGGGGMFGSGGGGGGTGSTGPGYSFPHYGFPTYGGITFHPGTTKSNAGMKHGTMDTESKKDPEGCDKSDDKNTVNLFGKVIETTEQDQEPSEATVGNGEVTLTYATGTKEESAGVQDNLFLEKAMQLAKRHANALFDYAVTGDVKMLLAVQRHLTAVQDENGDSVLHLAIIHLHSQLVRDLLEVTSGLISDDIINMRNDLYQTPLHLAVITKQEDVVEDLLRAGADLSLLDRLGNSVLHLAAKEGHDKVLSILLKHKKAALLLDHPNGDGLNAIHLAMMSNSLPCLLLLVAAGADVNAQEQKSGRTALHLAVEHDNISLAGCLLLEGDAHVDSTTYDGTTPLHIAAGRGSTRLAALLKAAGADPLVENFEPLYDLDDSWENAGEDEGVVPGTTPLDMATSWQVFDILNGKPYEPEFTSDDLLAQGDMKQLAEDVKLQLYKLLEIPDPDKNWATLAQKLGLGILNNAFRLSPAPSKTLMDNYEVSGGTVRELVEALRQMGYTEAIEVIQAASSPVKTTSQAHSLPLSPASTRQQIDELRDSDSVCDSGVETSFRKLSFTESLTSGASLLTLNKMPHDYGQEGPLEGKIXXMAEDDPYLGRPEQMFHLDPSLTHTIFNPEVFQPQMALPTADGPYLQILEQPKQRGFRFRYVCEGPSHGGLPGASSEKNKKSYPQVKICNYVGPAKVIVQLVTNGKNIHLHAHSLVGKHCEDGICTVTAGPKDMVVGFANLGILHVTKKKVFETLEARMTEACIRGYNPGLLVHPDLAYLQAEGGGDRQLGDREKELIRQAALQQTKEMDLSVVRLMFTAFLPDSTGSFTRRLEPVVSDAIYDSKAPNASNLKIVRMDRTAGCVTGGEEIYLLCDKVQKDDIQIRFYEEEENGGVWEGFGDFSPTDVHRQFAIVFKTPKYKDINITKPASVFVQLRRKSDLETSEPKPFLYYPEIKDKEEVQRKRQKLMPNFSDSFGGGSGAGAGGGGMFGSGGGGGGTGSTGPGYSFPHYGFPTYGGITFHPGTTKSNAGMKHGTMDTESKKDPEGCDKSDDKNTVNLFGKVIETTEQDQEPSEATVGNGEVTLTYATGTKEESAGVQDNLFLEKAMQLAKRHANALFDYAVTGDVKMLLAVQRHLTAVQDENGDSVLHLAIIHLHSQLVRDLLEVTSGLISDDIINMRNDLYQTPLHLAVITKQEDVVEDLLRAGADLSLLDRLGNSVLHLAAKEGHDKVLSILLKHKKAALLLDHPNGDGLNAIHLAMMSNSLPCLLLLVAAGADVNAQEQKSGRTALHLAVEHDNISLAGCLLLEGDAHVDSTTYDGTTPLHIAAGRGSTRLAALLKAAGADPLVENFEPLYDLDDSWENAGEDEGVVPGTTPLDMATSWQVFDILNGKPYEPEFTSDDLLAQGDMKQLAEDVKLQLYKLLEIPDPDKNWATLAQKLGLGILNNAFRLSPAPSKTLMDNYEVSGGTVRELVEALRQMGYTEAIEVIQAASSPVKTTSQAHSLPLSPASTRQQIDELRDSDSVCDSGVETSFRKLSFTESLTSGASLLTLNKMPHDYGQEGPLEGKIXXMAEDDPYLGRPEQMFHLDPSLTHTIFNPEVFQPQMALPTDGPYLQILEQPKQRGFRFRYVCEGPSHGGLPGASSEKNKKSYPQVKICNYVGPAKVIVQLVTNGKNIHLHAHSLVGKHCEDGICTVTAGPKDMVVGFANLGILHVTKKKVFETLEARMTEACIRGYNPGLLVHPDLAYLQAEGGGDRQLGDREKELIRQAALQQTKEMDLSVVRLMFTAFLPDSTGSFTRRLEPVVSDAIYDSKAPNASNLKIVRMDRTAGCVTGGEEIYLLCDKVQKDDIQIRFYEEEENGGVWEGFGDFSPTDVHRQFAIVFKTPKYKDINITKPASVFVQLRRKSDLETSEPKPFLYYPEIKDKEEVQRKRQKLMPNFSDSFGGGSGAGAGGGGMFGSGGGGGGTGSTGPGYSFPHYGFPTYGGITFHPGTTKSNAGMKHGTMDTESKKDPEGCDKSDDKNTVNLFGKVIETTEQDQEPSEATVGNGEVTLTYATGTKEESAGVQDNLFLEKAMQLAKRHANALFDYAVTGDVKMLLAVQRHLTAVQDENGDSVLHLAIIHLHSQLVRDLLEVTSGLISDDIINMRNDLYQTPLHLAVITKQEDVVEDLLRAGADLSLLDRLGNSVLHLAAKEGHDKVLSILLKHKKAALLLDHPNGDGLNAIHLAMMSNSLPCLLLLVAAGADVNAQEQKSGRTALHLAVEHDNISLAGCLLLEGDAHVDSTTYDGTTPLHIAAGRGSTRLAALLKAAGADPLVENFEPLYDLDDSWENAGEDEGVVPGTTPLDMATSWQVFDILNGKPYEPEFTSDDLLAQGDMKQLAEDVKLQLYKLLEIPDPDKNWATLAQKLGLGILNNAFRLSPAPSKTLMDNYEVSGGTVRELVEALRQMGYTEAIEVIQAASSPVKTTSQAHSLPLSPASTRQQIDELRDSDSVCDSGVETSFRKLSFTESLTSGASLLTLNKMPHDYGQEGPLEGKIXXMFHLDPSLTHTIFNPEVFQPQMALPTDGPYLQILEQPKQRGFRFRYVCEGPSHGGLPGASSEKNKKSYPQVKICNYVGPAKVIVQLVTNGKNIHLHAHSLVGKHCEDGICTVTAGPKDMVVGFANLGILHVTKKKVFETLEARMTEACIRGYNPGLLVHPDLAYLQAEGGGDRQLGDREKELIRQAALQQTKEMDLSVVRLMFTAFLPDSTGSFTRRLEPVVSDAIYDSKAPNASNLKIVRMDRTAGCVTGGEEIYLLCDKVQKDDIQIRFYEEEENGGVWEGFGDFSPTDVHRQFAIVFKTPKYKDINITKPASVFVQLRRKSDLETSEPKPFLYYPEIKDKEEVQRKRQKLMPNFSDSFGGGSGAGAGGGGMFGSGGGGGGTGSTGPGYSFPHYGFPTYGGITFHPGTTKSNAGMKHGTMDTESKKDPEGCDKSDDKNTVNLFGKVIETTEQDQEPSEATVGNGEVTLTYATGTKEESAGVQDNLFLEKAMQLAKRHANALFDYAVTGDVKMLLAVQRHLTAVQDENGDSVLHLAIIHLHSQLVRDLLEVTSGLISDDIINMRNDLYQTPLHLAVITKQEDVVEDLLRAGADLSLLDRLGNSVLHLAAKEGHDKVLSILLKHKKAALLLDHPNGDGLNAIHLAMMSNSLPCLLLLVAAGADVNAQEQKSGRTALHLAVEHDNISLAGCLLLEGDAHVDSTTYDGTTPLHIAAGRGSTRLAALLKAAGADPLVENFEPLYDLDDSWENAGEDEGVVPGTTPLDMATSWQVFDILNGKPYEPEFTSDDLLAQGDMKQLAEDVKLQLYKLLEIPDPDKNWATLAQKLGLGILNNAFRLSPAPSKTLMDNYEVSGGTVRELVEALRQMGYTEAIEVIQAASSPVKTTSQAHSLPLSPASTRQQIDELRDSDSVCDSGVETSFRKLSFTESLTSGASLLTLNKMPHDYGQEGPLEGKIXXMAEDDPYLGRPEQMFHLDPSLTHTIFNPEVFQPQMALPTADGPYLQILEQPKQRGFRFRYVCEGPSHGGLPGASSEKNKKSYPQVKICNYVGPAKVIVQLVTNGKNIHLHAHSLVGKHCEDGICTVTAGPKDMVVGFANLGILHVTKKKVFETLEARMTEACIRGYNPGLLVHPDLAYLQAEGGGDRQLGDREKELIRQAALQQTKEMDLSVVRLMFTAFLPDSTGSFTRRLEPVVSDAIYDSKAPNASNLKIVRMDRTAGCVTGGEEIYLLCDKVQKDDIQIRFYEEEENGGVWEGFGDFSPTDVHRQFAIVFKTPKYKDINITKPASVFVQLRRKSDLETSEPKPFLYYPEIKDKEEVQRKRQKLMPNFSDSFGGGSGAGAGGGGMFGSGGGGGGTGSTGPGYSFPHYGFPTYGGITFHPGTTKSNAGMKHGTMDTESKKDPEGCDKSDDKNTVNLFGKVIETTEQDQEPSEATVGNGEVTLTYATGTKEESAGVQDNLFLEKAMQLAKRHANALFDYAVTGDVKMLLAVQRHLTAVQDENGDSVLHLAIIHLHSQLVRDLLEVTSGLISDDIINMRNDLYQTPLHLAVITKQEDVVEDLLRAGADLSLLDRLGNSVLHLAAKEGHDKVLSILLKHKKAALLLDHPNGDGLNAIHLAMMSNSLPCLLLLVAAGADVNAQEQKSGRTALHLAVEHDNISLAGCLLLEGDAHVDSTTYDGTTPLHIAAGRGSTRLAALLKAAGADPLVENFEPLYDLDDSWENAGEDEGVVPGTTPLDMATSWQVFDILNGKPYEPEFTSDDLLAQGDMKQLAEDVKLQLYKLLEIPDPDKNWATLAQKLGLGILNNAFRLSPAPSKTLMDNYEVSGGTVRELVEALRQMGYTEAIEVIQAASSPVKTTSQAHSLPLSPASTRQQIDELRDSDSVCDSGVETSFRKLSFTESLTSGASLLTLNKMPHDYGQEGPLEGKIX"
## [3,] "MPRCPAGAMDEGPVDLRTRPKAAGLPGAALPLRKRPLRAPSPEPAAPRGAAGLVVPLDPLRGGCDLPAVPGPPHGLARPEALYYPGALLPLYPTRAMGSPFPLVNLPTPLYPMMCPMEHPLSADIAMATRADEDGDTPLHIAVVQGNLPAVHRLVNLFQQGGRELDIYNNLRQTPLHLAVITTLPSVVRLLVTAGASPMALDRHGQTAAHLACEHRSPTCLRALLDSAAPGTLDLEARNYDGLTALHVAVNTECQETVQLLLERGADIDAVDIKSGRSPLIHAVENNSLSMVQLLLQHGANVNAQMYSGSSALHSASGRGLLPLVRTLVRSGADSSLKNCHNDTPLMVARSRRVIDILRGKATRPASTSQPDPSPDRSANTSPESSSRLSSNGLLSASPSSSPSQSPPRDPPGFPMAPPNFFLPSPSPPAFLPFAGVLRGPGRPVPPSPAPGGSXXMAEDDPYLGRPEQMFHLDPSLTHTIFNPEVFQPQMALPTDGPYLQILEQPKQRGFRFRYVCEGPSHGGLPGASSEKNKKSYPQVKICNYVGPAKVIVQLVTNGKNIHLHAHSLVGKHCEDGICTVTAGPKDMVVGFANLGILHVTKKKVFETLEARMTEACIRGYNPGLLVHPDLAYLQAEGGGDRQLGDREKELIRQAALQQTKEMDLSVVRLMFTAFLPDSTGSFTRRLEPVVSDAIYDSKAPNASNLKIVRMDRTAGCVTGGEEIYLLCDKVQKDDIQIRFYEEEENGGVWEGFGDFSPTDVHRQFAIVFKTPKYKDINITKPASVFVQLRRKSDLETSEPKPFLYYPEIKDKEEVQRKRQKLMPNFSDSFGGGSGAGAGGGGMFGSGGGGGGTGSTGPGYSFPHYGFPTYGGITFHPGTTKSNAGMKHGTMDTESKKDPEGCDKSDDKNTVNLFGKVIETTEQDQEPSEATVGNGEVTLTYATGTKEESAGVQDNLFLEKAMQLAKRHANALFDYAVTGDVKMLLAVQRHLTAVQDENGDSVLHLAIIHLHSQLVRDLLEVTSGLISDDIINMRNDLYQTPLHLAVITKQEDVVEDLLRAGADLSLLDRLGNSVLHLAAKEGHDKVLSILLKHKKAALLLDHPNGDGLNAIHLAMMSNSLPCLLLLVAAGADVNAQEQKSGRTALHLAVEHDNISLAGCLLLEGDAHVDSTTYDGTTPLHIAAGRGSTRLAALLKAAGADPLVENFEPLYDLDDSWENAGEDEGVVPGTTPLDMATSWQVFDILNGKPYEPEFTSDDLLAQGDMKQLAEDVKLQLYKLLEIPDPDKNWATLAQKLGLGILNNAFRLSPAPSKTLMDNYEVSGGTVRELVEALRQMGYTEAIEVIQAASSPVKTTSQAHSLPLSPASTRQQIDELRDSDSVCDSGVETSFRKLSFTESLTSGASLLTLNKMPHDYGQEGPLEGKIXXMAEDDPYLGRPEQMFHLDPSLTHTIFNPEVFQPQMALPTADGPYLQILEQPKQRGFRFRYVCEGPSHGGLPGASSEKNKKSYPQVKICNYVGPAKVIVQLVTNGKNIHLHAHSLVGKHCEDGICTVTAGPKDMVVGFANLGILHVTKKKVFETLEARMTEACIRGYNPGLLVHPDLAYLQAEGGGDRQLGDREKELIRQAALQQTKEMDLSVVRLMFTAFLPDSTGSFTRRLEPVVSDAIYDSKAPNASNLKIVRMDRTAGCVTGGEEIYLLCDKVQKDDIQIRFYEEEENGGVWEGFGDFSPTDVHRQFAIVFKTPKYKDINITKPASVFVQLRRKSDLETSEPKPFLYYPEIKDKEEVQRKRQKLMPNFSDSFGGGSGAGAGGGGMFGSGGGGGGTGSTGPGYSFPHYGFPTYGGITFHPGTTKSNAGMKHGTMDTESKKDPEGCDKSDDKNTVNLFGKVIETTEQDQEPSEATVGNGEVTLTYATGTKEESAGVQDNLFLEKAMQLAKRHANALFDYAVTGDVKMLLAVQRHLTAVQDENGDSVLHLAIIHLHSQLVRDLLEVTSGLISDDIINMRNDLYQTPLHLAVITKQEDVVEDLLRAGADLSLLDRLGNSVLHLAAKEGHDKVLSILLKHKKAALLLDHPNGDGLNAIHLAMMSNSLPCLLLLVAAGADVNAQEQKSGRTALHLAVEHDNISLAGCLLLEGDAHVDSTTYDGTTPLHIAAGRGSTRLAALLKAAGADPLVENFEPLYDLDDSWENAGEDEGVVPGTTPLDMATSWQVFDILNGKPYEPEFTSDDLLAQGDMKQLAEDVKLQLYKLLEIPDPDKNWATLAQKLGLGILNNAFRLSPAPSKTLMDNYEVSGGTVRELVEALRQMGYTEAIEVIQAASSPVKTTSQAHSLPLSPASTRQQIDELRDSDSVCDSGVETSFRKLSFTESLTSGASLLTLNKMPHDYGQEGPLEGKIXXMAEDDPYLGRPEQMFHLDPSLTHTIFNPEVFQPQMALPTDGPYLQILEQPKQRGFRFRYVCEGPSHGGLPGASSEKNKKSYPQVKICNYVGPAKVIVQLVTNGKNIHLHAHSLVGKHCEDGICTVTAGPKDMVVGFANLGILHVTKKKVFETLEARMTEACIRGYNPGLLVHPDLAYLQAEGGGDRQLGDREKELIRQAALQQTKEMDLSVVRLMFTAFLPDSTGSFTRRLEPVVSDAIYDSKAPNASNLKIVRMDRTAGCVTGGEEIYLLCDKVQKDDIQIRFYEEEENGGVWEGFGDFSPTDVHRQFAIVFKTPKYKDINITKPASVFVQLRRKSDLETSEPKPFLYYPEIKDKEEVQRKRQKLMPNFSDSFGGGSGAGAGGGGMFGSGGGGGGTGSTGPGYSFPHYGFPTYGGITFHPGTTKSNAGMKHGTMDTESKKDPEGCDKSDDKNTVNLFGKVIETTEQDQEPSEATVGNGEVTLTYATGTKEESAGVQDNLFLEKAMQLAKRHANALFDYAVTGDVKMLLAVQRHLTAVQDENGDSVLHLAIIHLHSQLVRDLLEVTSGLISDDIINMRNDLYQTPLHLAVITKQEDVVEDLLRAGADLSLLDRLGNSVLHLAAKEGHDKVLSILLKHKKAALLLDHPNGDGLNAIHLAMMSNSLPCLLLLVAAGADVNAQEQKSGRTALHLAVEHDNISLAGCLLLEGDAHVDSTTYDGTTPLHIAAGRGSTRLAALLKAAGADPLVENFEPLYDLDDSWENAGEDEGVVPGTTPLDMATSWQVFDILNGKPYEPEFTSDDLLAQGDMKQLAEDVKLQLYKLLEIPDPDKNWATLAQKLGLGILNNAFRLSPAPSKTLMDNYEVSGGTVRELVEALRQMGYTEAIEVIQAASSPVKTTSQAHSLPLSPASTRQQIDELRDSDSVCDSGVETSFRKLSFTESLTSGASLLTLNKMPHDYGQEGPLEGKIXXMAEDDPYLGRPEQMFHLDPSLTHTIFNPEVFQPQMALPTADGPYLQILEQPKQRGFRFRYVCEGPSHGGLPGASSEKNKKSYPQVKICNYVGPAKVIVQLVTNGKNIHLHAHSLVGKHCEDGICTVTAGPKDMVVGFANLGILHVTKKKVFETLEARMTEACIRGYNPGLLVHPDLAYLQAEGGGDRQLGDREKELIRQAALQQTKEMDLSVVRLMFTAFLPDSTGSFTRRLEPVVSDAIYDSKAPNASNLKIVRMDRTAGCVTGGEEIYLLCDKVQKDDIQIRFYEEEENGGVWEGFGDFSPTDVHRQFAIVFKTPKYKDINITKPASVFVQLRRKSDLETSEPKPFLYYPEIKDKEEVQRKRQKLMPNFSDSFGGGSGAGAGGGGMFGSGGGGGGTGSTGPGYSFPHYGFPTYGGITFHPGTTKSNAGMKHGTMDTESKKDPEGCDKSDDKNTVNLFGKVIETTEQDQEPSEATVGNGEVTLTYATGTKEESAGVQDNLFLEKAMQLAKRHANALFDYAVTGDVKMLLAVQRHLTAVQDENGDSVLHLAIIHLHSQLVRDLLEVTSGLISDDIINMRNDLYQTPLHLAVITKQEDVVEDLLRAGADLSLLDRLGNSVLHLAAKEGHDKVLSILLKHKKAALLLDHPNGDGLNAIHLAMMSNSLPCLLLLVAAGADVNAQEQKSGRTALHLAVEHDNISLAGCLLLEGDAHVDSTTYDGTTPLHIAAGRGSTRLAALLKAAGADPLVENFEPLYDLDDSWENAGEDEGVVPGTTPLDMATSWQVFDILNGKPYEPEFTSDDLLAQGDMKQLAEDVKLQLYKLLEIPDPDKNWATLAQKLGLGILNNAFRLSPAPSKTLMDNYEVSGGTVRELVEALRQMGYTEAIEVIQAASSPVKTTSQAHSLPLSPASTRQQIDELRDSDSVCDSGVETSFRKLSFTESLTSGASLLTLNKMPHDYGQEGPLEGKIXXMAEDDPYLGRPEQMFHLDPSLTHTIFNPEVFQPQMALPTDGPYLQILEQPKQRGFRFRYVCEGPSHGGLPGASSEKNKKSYPQVKICNYVGPAKVIVQLVTNGKNIHLHAHSLVGKHCEDGICTVTAGPKDMVVGFANLGILHVTKKKVFETLEARMTEACIRGYNPGLLVHPDLAYLQAEGGGDRQLGDREKELIRQAALQQTKEMDLSVVRLMFTAFLPDSTGSFTRRLEPVVSDAIYDSKAPNASNLKIVRMDRTAGCVTGGEEIYLLCDKVQKDDIQIRFYEEEENGGVWEGFGDFSPTDVHRQFAIVFKTPKYKDINITKPASVFVQLRRKSDLETSEPKPFLYYPEIKDKEEVQRKRQKLMPNFSDSFGGGSGAGAGGGGMFGSGGGGGGTGSTGPGYSFPHYGFPTYGGITFHPGTTKSNAGMKHGTMDTESKKDPEGCDKSDDKNTVNLFGKVIETTEQDQEPSEATVGNGEVTLTYATGTKEESAGVQDNLFLEKAMQLAKRHANALFDYAVTGDVKMLLAVQRHLTAVQDENGDSVLHLAIIHLHSQLVRDLLEVTSGLISDDIINMRNDLYQTPLHLAVITKQEDVVEDLLRAGADLSLLDRLGNSVLHLAAKEGHDKVLSILLKHKKAALLLDHPNGDGLNAIHLAMMSNSLPCLLLLVAAGADVNAQEQKSGRTALHLAVEHDNISLAGCLLLEGDAHVDSTTYDGTTPLHIAAGRGSTRLAALLKAAGADPLVENFEPLYDLDDSWENAGEDEGVVPGTTPLDMATSWQVFDILNGKPYEPEFTSDDLLAQGDMKQLAEDVKLQLYKLLEIPDPDKNWATLAQKLGLGILNNAFRLSPAPSKTLMDNYEVSGGTVRELVEALRQMGYTEAIEVIQAASSPVKTTSQAHSLPLSPASTRQQIDELRDSDSVCDSGVETSFRKLSFTESLTSGASLLTLNKMPHDYGQEGPLEGKIXXMFHLDPSLTHTIFNPEVFQPQMALPTDGPYLQILEQPKQRGFRFRYVCEGPSHGGLPGASSEKNKKSYPQVKICNYVGPAKVIVQLVTNGKNIHLHAHSLVGKHCEDGICTVTAGPKDMVVGFANLGILHVTKKKVFETLEARMTEACIRGYNPGLLVHPDLAYLQAEGGGDRQLGDREKELIRQAALQQTKEMDLSVVRLMFTAFLPDSTGSFTRRLEPVVSDAIYDSKAPNASNLKIVRMDRTAGCVTGGEEIYLLCDKVQKDDIQIRFYEEEENGGVWEGFGDFSPTDVHRQFAIVFKTPKYKDINITKPASVFVQLRRKSDLETSEPKPFLYYPEIKDKEEVQRKRQKLMPNFSDSFGGGSGAGAGGGGMFGSGGGGGGTGSTGPGYSFPHYGFPTYGGITFHPGTTKSNAGMKHGTMDTESKKDPEGCDKSDDKNTVNLFGKVIETTEQDQEPSEATVGNGEVTLTYATGTKEESAGVQDNLFLEKAMQLAKRHANALFDYAVTGDVKMLLAVQRHLTAVQDENGDSVLHLAIIHLHSQLVRDLLEVTSGLISDDIINMRNDLYQTPLHLAVITKQEDVVEDLLRAGADLSLLDRLGNSVLHLAAKEGHDKVLSILLKHKKAALLLDHPNGDGLNAIHLAMMSNSLPCLLLLVAAGADVNAQEQKSGRTALHLAVEHDNISLAGCLLLEGDAHVDSTTYDGTTPLHIAAGRGSTRLAALLKAAGADPLVENFEPLYDLDDSWENAGEDEGVVPGTTPLDMATSWQVFDILNGKPYEPEFTSDDLLAQGDMKQLAEDVKLQLYKLLEIPDPDKNWATLAQKLGLGILNNAFRLSPAPSKTLMDNYEVSGGTVRELVEALRQMGYTEAIEVIQAASSPVKTTSQAHSLPLSPASTRQQIDELRDSDSVCDSGVETSFRKLSFTESLTSGASLLTLNKMPHDYGQEGPLEGKIXXMAEDDPYLGRPEQMFHLDPSLTHTIFNPEVFQPQMALPTADGPYLQILEQPKQRGFRFRYVCEGPSHGGLPGASSEKNKKSYPQVKICNYVGPAKVIVQLVTNGKNIHLHAHSLVGKHCEDGICTVTAGPKDMVVGFANLGILHVTKKKVFETLEARMTEACIRGYNPGLLVHPDLAYLQAEGGGDRQLGDREKELIRQAALQQTKEMDLSVVRLMFTAFLPDSTGSFTRRLEPVVSDAIYDSKAPNASNLKIVRMDRTAGCVTGGEEIYLLCDKVQKDDIQIRFYEEEENGGVWEGFGDFSPTDVHRQFAIVFKTPKYKDINITKPASVFVQLRRKSDLETSEPKPFLYYPEIKDKEEVQRKRQKLMPNFSDSFGGGSGAGAGGGGMFGSGGGGGGTGSTGPGYSFPHYGFPTYGGITFHPGTTKSNAGMKHGTMDTESKKDPEGCDKSDDKNTVNLFGKVIETTEQDQEPSEATVGNGEVTLTYATGTKEESAGVQDNLFLEKAMQLAKRHANALFDYAVTGDVKMLLAVQRHLTAVQDENGDSVLHLAIIHLHSQLVRDLLEVTSGLISDDIINMRNDLYQTPLHLAVITKQEDVVEDLLRAGADLSLLDRLGNSVLHLAAKEGHDKVLSILLKHKKAALLLDHPNGDGLNAIHLAMMSNSLPCLLLLVAAGADVNAQEQKSGRTALHLAVEHDNISLAGCLLLEGDAHVDSTTYDGTTPLHIAAGRGSTRLAALLKAAGADPLVENFEPLYDLDDSWENAGEDEGVVPGTTPLDMATSWQVFDILNGKPYEPEFTSDDLLAQGDMKQLAEDVKLQLYKLLEIPDPDKNWATLAQKLGLGILNNAFRLSPAPSKTLMDNYEVSGGTVRELVEALRQMGYTEAIEVIQAASSPVKTTSQAHSLPLSPASTRQQIDELRDSDSVCDSGVETSFRKLSFTESLTSGASLLTLNKMPHDYGQEGPLEGKIX"
## [4,] "MPRCPAGAMDEGPVDLRTRPKAAGLPGAALPLRKRPLRAPSPEPAAPRGAAGLVVPLDPLRGGCDLPAVPGPPHGLARPEALYYPGALLPLYPTRAMGSPFPLVNLPTPLYPMMCPMEHPLSADIAMATRADEDGDTPLHIAVVQGNLPAVHRLVNLFQQGGRELDIYNNLRQTPLHLAVITTLPSVVRLLVTAGASPMALDRHGQTAAHLACEHRSPTCLRALLDSAAPGTLDLEARNYDGLTALHVAVNTECQETVQLLLERGADIDAVDIKSGRSPLIHAVENNSLSMVQLLLQHGANVNAQMYSGSSALHSASGRGLLPLVRTLVRSGADSSLKNCHNDTPLMVARSRRVIDILRGKATRPASTSQPDPSPDRSANTSPESSSRLSSNGLLSASPSSSPSQSPPRDPPGFPMAPPNFFLPSPSPPAFLPFAGVLRGPGRPVPPSPAPGGSXXMAEDDPYLGRPEQMFHLDPSLTHTIFNPEVFQPQMALPTDGPYLQILEQPKQRGFRFRYVCEGPSHGGLPGASSEKNKKSYPQVKICNYVGPAKVIVQLVTNGKNIHLHAHSLVGKHCEDGICTVTAGPKDMVVGFANLGILHVTKKKVFETLEARMTEACIRGYNPGLLVHPDLAYLQAEGGGDRQLGDREKELIRQAALQQTKEMDLSVVRLMFTAFLPDSTGSFTRRLEPVVSDAIYDSKAPNASNLKIVRMDRTAGCVTGGEEIYLLCDKVQKDDIQIRFYEEEENGGVWEGFGDFSPTDVHRQFAIVFKTPKYKDINITKPASVFVQLRRKSDLETSEPKPFLYYPEIKDKEEVQRKRQKLMPNFSDSFGGGSGAGAGGGGMFGSGGGGGGTGSTGPGYSFPHYGFPTYGGITFHPGTTKSNAGMKHGTMDTESKKDPEGCDKSDDKNTVNLFGKVIETTEQDQEPSEATVGNGEVTLTYATGTKEESAGVQDNLFLEKAMQLAKRHANALFDYAVTGDVKMLLAVQRHLTAVQDENGDSVLHLAIIHLHSQLVRDLLEVTSGLISDDIINMRNDLYQTPLHLAVITKQEDVVEDLLRAGADLSLLDRLGNSVLHLAAKEGHDKVLSILLKHKKAALLLDHPNGDGLNAIHLAMMSNSLPCLLLLVAAGADVNAQEQKSGRTALHLAVEHDNISLAGCLLLEGDAHVDSTTYDGTTPLHIAAGRGSTRLAALLKAAGADPLVENFEPLYDLDDSWENAGEDEGVVPGTTPLDMATSWQVFDILNGKPYEPEFTSDDLLAQGDMKQLAEDVKLQLYKLLEIPDPDKNWATLAQKLGLGILNNAFRLSPAPSKTLMDNYEVSGGTVRELVEALRQMGYTEAIEVIQAASSPVKTTSQAHSLPLSPASTRQQIDELRDSDSVCDSGVETSFRKLSFTESLTSGASLLTLNKMPHDYGQEGPLEGKIXXMAEDDPYLGRPEQMFHLDPSLTHTIFNPEVFQPQMALPTADGPYLQILEQPKQRGFRFRYVCEGPSHGGLPGASSEKNKKSYPQVKICNYVGPAKVIVQLVTNGKNIHLHAHSLVGKHCEDGICTVTAGPKDMVVGFANLGILHVTKKKVFETLEARMTEACIRGYNPGLLVHPDLAYLQAEGGGDRQLGDREKELIRQAALQQTKEMDLSVVRLMFTAFLPDSTGSFTRRLEPVVSDAIYDSKAPNASNLKIVRMDRTAGCVTGGEEIYLLCDKVQKDDIQIRFYEEEENGGVWEGFGDFSPTDVHRQFAIVFKTPKYKDINITKPASVFVQLRRKSDLETSEPKPFLYYPEIKDKEEVQRKRQKLMPNFSDSFGGGSGAGAGGGGMFGSGGGGGGTGSTGPGYSFPHYGFPTYGGITFHPGTTKSNAGMKHGTMDTESKKDPEGCDKSDDKNTVNLFGKVIETTEQDQEPSEATVGNGEVTLTYATGTKEESAGVQDNLFLEKAMQLAKRHANALFDYAVTGDVKMLLAVQRHLTAVQDENGDSVLHLAIIHLHSQLVRDLLEVTSGLISDDIINMRNDLYQTPLHLAVITKQEDVVEDLLRAGADLSLLDRLGNSVLHLAAKEGHDKVLSILLKHKKAALLLDHPNGDGLNAIHLAMMSNSLPCLLLLVAAGADVNAQEQKSGRTALHLAVEHDNISLAGCLLLEGDAHVDSTTYDGTTPLHIAAGRGSTRLAALLKAAGADPLVENFEPLYDLDDSWENAGEDEGVVPGTTPLDMATSWQVFDILNGKPYEPEFTSDDLLAQGDMKQLAEDVKLQLYKLLEIPDPDKNWATLAQKLGLGILNNAFRLSPAPSKTLMDNYEVSGGTVRELVEALRQMGYTEAIEVIQAASSPVKTTSQAHSLPLSPASTRQQIDELRDSDSVCDSGVETSFRKLSFTESLTSGASLLTLNKMPHDYGQEGPLEGKIXXMAEDDPYLGRPEQMFHLDPSLTHTIFNPEVFQPQMALPTDGPYLQILEQPKQRGFRFRYVCEGPSHGGLPGASSEKNKKSYPQVKICNYVGPAKVIVQLVTNGKNIHLHAHSLVGKHCEDGICTVTAGPKDMVVGFANLGILHVTKKKVFETLEARMTEACIRGYNPGLLVHPDLAYLQAEGGGDRQLGDREKELIRQAALQQTKEMDLSVVRLMFTAFLPDSTGSFTRRLEPVVSDAIYDSKAPNASNLKIVRMDRTAGCVTGGEEIYLLCDKVQKDDIQIRFYEEEENGGVWEGFGDFSPTDVHRQFAIVFKTPKYKDINITKPASVFVQLRRKSDLETSEPKPFLYYPEIKDKEEVQRKRQKLMPNFSDSFGGGSGAGAGGGGMFGSGGGGGGTGSTGPGYSFPHYGFPTYGGITFHPGTTKSNAGMKHGTMDTESKKDPEGCDKSDDKNTVNLFGKVIETTEQDQEPSEATVGNGEVTLTYATGTKEESAGVQDNLFLEKAMQLAKRHANALFDYAVTGDVKMLLAVQRHLTAVQDENGDSVLHLAIIHLHSQLVRDLLEVTSGLISDDIINMRNDLYQTPLHLAVITKQEDVVEDLLRAGADLSLLDRLGNSVLHLAAKEGHDKVLSILLKHKKAALLLDHPNGDGLNAIHLAMMSNSLPCLLLLVAAGADVNAQEQKSGRTALHLAVEHDNISLAGCLLLEGDAHVDSTTYDGTTPLHIAAGRGSTRLAALLKAAGADPLVENFEPLYDLDDSWENAGEDEGVVPGTTPLDMATSWQVFDILNGKPYEPEFTSDDLLAQGDMKQLAEDVKLQLYKLLEIPDPDKNWATLAQKLGLGILNNAFRLSPAPSKTLMDNYEVSGGTVRELVEALRQMGYTEAIEVIQAASSPVKTTSQAHSLPLSPASTRQQIDELRDSDSVCDSGVETSFRKLSFTESLTSGASLLTLNKMPHDYGQEGPLEGKIXXMAEDDPYLGRPEQMFHLDPSLTHTIFNPEVFQPQMALPTADGPYLQILEQPKQRGFRFRYVCEGPSHGGLPGASSEKNKKSYPQVKICNYVGPAKVIVQLVTNGKNIHLHAHSLVGKHCEDGICTVTAGPKDMVVGFANLGILHVTKKKVFETLEARMTEACIRGYNPGLLVHPDLAYLQAEGGGDRQLGDREKELIRQAALQQTKEMDLSVVRLMFTAFLPDSTGSFTRRLEPVVSDAIYDSKAPNASNLKIVRMDRTAGCVTGGEEIYLLCDKVQKDDIQIRFYEEEENGGVWEGFGDFSPTDVHRQFAIVFKTPKYKDINITKPASVFVQLRRKSDLETSEPKPFLYYPEIKDKEEVQRKRQKLMPNFSDSFGGGSGAGAGGGGMFGSGGGGGGTGSTGPGYSFPHYGFPTYGGITFHPGTTKSNAGMKHGTMDTESKKDPEGCDKSDDKNTVNLFGKVIETTEQDQEPSEATVGNGEVTLTYATGTKEESAGVQDNLFLEKAMQLAKRHANALFDYAVTGDVKMLLAVQRHLTAVQDENGDSVLHLAIIHLHSQLVRDLLEVTSGLISDDIINMRNDLYQTPLHLAVITKQEDVVEDLLRAGADLSLLDRLGNSVLHLAAKEGHDKVLSILLKHKKAALLLDHPNGDGLNAIHLAMMSNSLPCLLLLVAAGADVNAQEQKSGRTALHLAVEHDNISLAGCLLLEGDAHVDSTTYDGTTPLHIAAGRGSTRLAALLKAAGADPLVENFEPLYDLDDSWENAGEDEGVVPGTTPLDMATSWQVFDILNGKPYEPEFTSDDLLAQGDMKQLAEDVKLQLYKLLEIPDPDKNWATLAQKLGLGILNNAFRLSPAPSKTLMDNYEVSGGTVRELVEALRQMGYTEAIEVIQAASSPVKTTSQAHSLPLSPASTRQQIDELRDSDSVCDSGVETSFRKLSFTESLTSGASLLTLNKMPHDYGQEGPLEGKIXXMAEDDPYLGRPEQMFHLDPSLTHTIFNPEVFQPQMALPTDGPYLQILEQPKQRGFRFRYVCEGPSHGGLPGASSEKNKKSYPQVKICNYVGPAKVIVQLVTNGKNIHLHAHSLVGKHCEDGICTVTAGPKDMVVGFANLGILHVTKKKVFETLEARMTEACIRGYNPGLLVHPDLAYLQAEGGGDRQLGDREKELIRQAALQQTKEMDLSVVRLMFTAFLPDSTGSFTRRLEPVVSDAIYDSKAPNASNLKIVRMDRTAGCVTGGEEIYLLCDKVQKDDIQIRFYEEEENGGVWEGFGDFSPTDVHRQFAIVFKTPKYKDINITKPASVFVQLRRKSDLETSEPKPFLYYPEIKDKEEVQRKRQKLMPNFSDSFGGGSGAGAGGGGMFGSGGGGGGTGSTGPGYSFPHYGFPTYGGITFHPGTTKSNAGMKHGTMDTESKKDPEGCDKSDDKNTVNLFGKVIETTEQDQEPSEATVGNGEVTLTYATGTKEESAGVQDNLFLEKAMQLAKRHANALFDYAVTGDVKMLLAVQRHLTAVQDENGDSVLHLAIIHLHSQLVRDLLEVTSGLISDDIINMRNDLYQTPLHLAVITKQEDVVEDLLRAGADLSLLDRLGNSVLHLAAKEGHDKVLSILLKHKKAALLLDHPNGDGLNAIHLAMMSNSLPCLLLLVAAGADVNAQEQKSGRTALHLAVEHDNISLAGCLLLEGDAHVDSTTYDGTTPLHIAAGRGSTRLAALLKAAGADPLVENFEPLYDLDDSWENAGEDEGVVPGTTPLDMATSWQVFDILNGKPYEPEFTSDDLLAQGDMKQLAEDVKLQLYKLLEIPDPDKNWATLAQKLGLGILNNAFRLSPAPSKTLMDNYEVSGGTVRELVEALRQMGYTEAIEVIQAASSPVKTTSQAHSLPLSPASTRQQIDELRDSDSVCDSGVETSFRKLSFTESLTSGASLLTLNKMPHDYGQEGPLEGKIXXMFHLDPSLTHTIFNPEVFQPQMALPTDGPYLQILEQPKQRGFRFRYVCEGPSHGGLPGASSEKNKKSYPQVKICNYVGPAKVIVQLVTNGKNIHLHAHSLVGKHCEDGICTVTAGPKDMVVGFANLGILHVTKKKVFETLEARMTEACIRGYNPGLLVHPDLAYLQAEGGGDRQLGDREKELIRQAALQQTKEMDLSVVRLMFTAFLPDSTGSFTRRLEPVVSDAIYDSKAPNASNLKIVRMDRTAGCVTGGEEIYLLCDKVQKDDIQIRFYEEEENGGVWEGFGDFSPTDVHRQFAIVFKTPKYKDINITKPASVFVQLRRKSDLETSEPKPFLYYPEIKDKEEVQRKRQKLMPNFSDSFGGGSGAGAGGGGMFGSGGGGGGTGSTGPGYSFPHYGFPTYGGITFHPGTTKSNAGMKHGTMDTESKKDPEGCDKSDDKNTVNLFGKVIETTEQDQEPSEATVGNGEVTLTYATGTKEESAGVQDNLFLEKAMQLAKRHANALFDYAVTGDVKMLLAVQRHLTAVQDENGDSVLHLAIIHLHSQLVRDLLEVTSGLISDDIINMRNDLYQTPLHLAVITKQEDVVEDLLRAGADLSLLDRLGNSVLHLAAKEGHDKVLSILLKHKKAALLLDHPNGDGLNAIHLAMMSNSLPCLLLLVAAGADVNAQEQKSGRTALHLAVEHDNISLAGCLLLEGDAHVDSTTYDGTTPLHIAAGRGSTRLAALLKAAGADPLVENFEPLYDLDDSWENAGEDEGVVPGTTPLDMATSWQVFDILNGKPYEPEFTSDDLLAQGDMKQLAEDVKLQLYKLLEIPDPDKNWATLAQKLGLGILNNAFRLSPAPSKTLMDNYEVSGGTVRELVEALRQMGYTEAIEVIQAASSPVKTTSQAHSLPLSPASTRQQIDELRDSDSVCDSGVETSFRKLSFTESLTSGASLLTLNKMPHDYGQEGPLEGKIXXMAEDDPYLGRPEQMFHLDPSLTHTIFNPEVFQPQMALPTADGPYLQILEQPKQRGFRFRYVCEGPSHGGLPGASSEKNKKSYPQVKICNYVGPAKVIVQLVTNGKNIHLHAHSLVGKHCEDGICTVTAGPKDMVVGFANLGILHVTKKKVFETLEARMTEACIRGYNPGLLVHPDLAYLQAEGGGDRQLGDREKELIRQAALQQTKEMDLSVVRLMFTAFLPDSTGSFTRRLEPVVSDAIYDSKAPNASNLKIVRMDRTAGCVTGGEEIYLLCDKVQKDDIQIRFYEEEENGGVWEGFGDFSPTDVHRQFAIVFKTPKYKDINITKPASVFVQLRRKSDLETSEPKPFLYYPEIKDKEEVQRKRQKLMPNFSDSFGGGSGAGAGGGGMFGSGGGGGGTGSTGPGYSFPHYGFPTYGGITFHPGTTKSNAGMKHGTMDTESKKDPEGCDKSDDKNTVNLFGKVIETTEQDQEPSEATVGNGEVTLTYATGTKEESAGVQDNLFLEKAMQLAKRHANALFDYAVTGDVKMLLAVQRHLTAVQDENGDSVLHLAIIHLHSQLVRDLLEVTSGLISDDIINMRNDLYQTPLHLAVITKQEDVVEDLLRAGADLSLLDRLGNSVLHLAAKEGHDKVLSILLKHKKAALLLDHPNGDGLNAIHLAMMSNSLPCLLLLVAAGADVNAQEQKSGRTALHLAVEHDNISLAGCLLLEGDAHVDSTTYDGTTPLHIAAGRGSTRLAALLKAAGADPLVENFEPLYDLDDSWENAGEDEGVVPGTTPLDMATSWQVFDILNGKPYEPEFTSDDLLAQGDMKQLAEDVKLQLYKLLEIPDPDKNWATLAQKLGLGILNNAFRLSPAPSKTLMDNYEVSGGTVRELVEALRQMGYTEAIEVIQAASSPVKTTSQAHSLPLSPASTRQQIDELRDSDSVCDSGVETSFRKLSFTESLTSGASLLTLNKMPHDYGQEGPLEGKIX"
## [5,] "MPRCPAGAMDEGPVDLRTRPKAAGLPGAALPLRKRPLRAPSPEPAAPRGAAGLVVPLDPLRGGCDLPAVPGPPHGLARPEALYYPGALLPLYPTRAMGSPFPLVNLPTPLYPMMCPMEHPLSADIAMATRADEDGDTPLHIAVVQGNLPAVHRLVNLFQQGGRELDIYNNLRQTPLHLAVITTLPSVVRLLVTAGASPMALDRHGQTAAHLACEHRSPTCLRALLDSAAPGTLDLEARNYDGLTALHVAVNTECQETVQLLLERGADIDAVDIKSGRSPLIHAVENNSLSMVQLLLQHGANVNAQMYSGSSALHSASGRGLLPLVRTLVRSGADSSLKNCHNDTPLMVARSRRVIDILRGKATRPASTSQPDPSPDRSANTSPESSSRLSSNGLLSASPSSSPSQSPPRDPPGFPMAPPNFFLPSPSPPAFLPFAGVLRGPGRPVPPSPAPGGSXXMAEDDPYLGRPEQMFHLDPSLTHTIFNPEVFQPQMALPTDGPYLQILEQPKQRGFRFRYVCEGPSHGGLPGASSEKNKKSYPQVKICNYVGPAKVIVQLVTNGKNIHLHAHSLVGKHCEDGICTVTAGPKDMVVGFANLGILHVTKKKVFETLEARMTEACIRGYNPGLLVHPDLAYLQAEGGGDRQLGDREKELIRQAALQQTKEMDLSVVRLMFTAFLPDSTGSFTRRLEPVVSDAIYDSKAPNASNLKIVRMDRTAGCVTGGEEIYLLCDKVQKDDIQIRFYEEEENGGVWEGFGDFSPTDVHRQFAIVFKTPKYKDINITKPASVFVQLRRKSDLETSEPKPFLYYPEIKDKEEVQRKRQKLMPNFSDSFGGGSGAGAGGGGMFGSGGGGGGTGSTGPGYSFPHYGFPTYGGITFHPGTTKSNAGMKHGTMDTESKKDPEGCDKSDDKNTVNLFGKVIETTEQDQEPSEATVGNGEVTLTYATGTKEESAGVQDNLFLEKAMQLAKRHANALFDYAVTGDVKMLLAVQRHLTAVQDENGDSVLHLAIIHLHSQLVRDLLEVTSGLISDDIINMRNDLYQTPLHLAVITKQEDVVEDLLRAGADLSLLDRLGNSVLHLAAKEGHDKVLSILLKHKKAALLLDHPNGDGLNAIHLAMMSNSLPCLLLLVAAGADVNAQEQKSGRTALHLAVEHDNISLAGCLLLEGDAHVDSTTYDGTTPLHIAAGRGSTRLAALLKAAGADPLVENFEPLYDLDDSWENAGEDEGVVPGTTPLDMATSWQVFDILNGKPYEPEFTSDDLLAQGDMKQLAEDVKLQLYKLLEIPDPDKNWATLAQKLGLGILNNAFRLSPAPSKTLMDNYEVSGGTVRELVEALRQMGYTEAIEVIQAASSPVKTTSQAHSLPLSPASTRQQIDELRDSDSVCDSGVETSFRKLSFTESLTSGASLLTLNKMPHDYGQEGPLEGKIXXMAEDDPYLGRPEQMFHLDPSLTHTIFNPEVFQPQMALPTADGPYLQILEQPKQRGFRFRYVCEGPSHGGLPGASSEKNKKSYPQVKICNYVGPAKVIVQLVTNGKNIHLHAHSLVGKHCEDGICTVTAGPKDMVVGFANLGILHVTKKKVFETLEARMTEACIRGYNPGLLVHPDLAYLQAEGGGDRQLGDREKELIRQAALQQTKEMDLSVVRLMFTAFLPDSTGSFTRRLEPVVSDAIYDSKAPNASNLKIVRMDRTAGCVTGGEEIYLLCDKVQKDDIQIRFYEEEENGGVWEGFGDFSPTDVHRQFAIVFKTPKYKDINITKPASVFVQLRRKSDLETSEPKPFLYYPEIKDKEEVQRKRQKLMPNFSDSFGGGSGAGAGGGGMFGSGGGGGGTGSTGPGYSFPHYGFPTYGGITFHPGTTKSNAGMKHGTMDTESKKDPEGCDKSDDKNTVNLFGKVIETTEQDQEPSEATVGNGEVTLTYATGTKEESAGVQDNLFLEKAMQLAKRHANALFDYAVTGDVKMLLAVQRHLTAVQDENGDSVLHLAIIHLHSQLVRDLLEVTSGLISDDIINMRNDLYQTPLHLAVITKQEDVVEDLLRAGADLSLLDRLGNSVLHLAAKEGHDKVLSILLKHKKAALLLDHPNGDGLNAIHLAMMSNSLPCLLLLVAAGADVNAQEQKSGRTALHLAVEHDNISLAGCLLLEGDAHVDSTTYDGTTPLHIAAGRGSTRLAALLKAAGADPLVENFEPLYDLDDSWENAGEDEGVVPGTTPLDMATSWQVFDILNGKPYEPEFTSDDLLAQGDMKQLAEDVKLQLYKLLEIPDPDKNWATLAQKLGLGILNNAFRLSPAPSKTLMDNYEVSGGTVRELVEALRQMGYTEAIEVIQAASSPVKTTSQAHSLPLSPASTRQQIDELRDSDSVCDSGVETSFRKLSFTESLTSGASLLTLNKMPHDYGQEGPLEGKIXXMAEDDPYLGRPEQMFHLDPSLTHTIFNPEVFQPQMALPTDGPYLQILEQPKQRGFRFRYVCEGPSHGGLPGASSEKNKKSYPQVKICNYVGPAKVIVQLVTNGKNIHLHAHSLVGKHCEDGICTVTAGPKDMVVGFANLGILHVTKKKVFETLEARMTEACIRGYNPGLLVHPDLAYLQAEGGGDRQLGDREKELIRQAALQQTKEMDLSVVRLMFTAFLPDSTGSFTRRLEPVVSDAIYDSKAPNASNLKIVRMDRTAGCVTGGEEIYLLCDKVQKDDIQIRFYEEEENGGVWEGFGDFSPTDVHRQFAIVFKTPKYKDINITKPASVFVQLRRKSDLETSEPKPFLYYPEIKDKEEVQRKRQKLMPNFSDSFGGGSGAGAGGGGMFGSGGGGGGTGSTGPGYSFPHYGFPTYGGITFHPGTTKSNAGMKHGTMDTESKKDPEGCDKSDDKNTVNLFGKVIETTEQDQEPSEATVGNGEVTLTYATGTKEESAGVQDNLFLEKAMQLAKRHANALFDYAVTGDVKMLLAVQRHLTAVQDENGDSVLHLAIIHLHSQLVRDLLEVTSGLISDDIINMRNDLYQTPLHLAVITKQEDVVEDLLRAGADLSLLDRLGNSVLHLAAKEGHDKVLSILLKHKKAALLLDHPNGDGLNAIHLAMMSNSLPCLLLLVAAGADVNAQEQKSGRTALHLAVEHDNISLAGCLLLEGDAHVDSTTYDGTTPLHIAAGRGSTRLAALLKAAGADPLVENFEPLYDLDDSWENAGEDEGVVPGTTPLDMATSWQVFDILNGKPYEPEFTSDDLLAQGDMKQLAEDVKLQLYKLLEIPDPDKNWATLAQKLGLGILNNAFRLSPAPSKTLMDNYEVSGGTVRELVEALRQMGYTEAIEVIQAASSPVKTTSQAHSLPLSPASTRQQIDELRDSDSVCDSGVETSFRKLSFTESLTSGASLLTLNKMPHDYGQEGPLEGKIXXMAEDDPYLGRPEQMFHLDPSLTHTIFNPEVFQPQMALPTADGPYLQILEQPKQRGFRFRYVCEGPSHGGLPGASSEKNKKSYPQVKICNYVGPAKVIVQLVTNGKNIHLHAHSLVGKHCEDGICTVTAGPKDMVVGFANLGILHVTKKKVFETLEARMTEACIRGYNPGLLVHPDLAYLQAEGGGDRQLGDREKELIRQAALQQTKEMDLSVVRLMFTAFLPDSTGSFTRRLEPVVSDAIYDSKAPNASNLKIVRMDRTAGCVTGGEEIYLLCDKVQKDDIQIRFYEEEENGGVWEGFGDFSPTDVHRQFAIVFKTPKYKDINITKPASVFVQLRRKSDLETSEPKPFLYYPEIKDKEEVQRKRQKLMPNFSDSFGGGSGAGAGGGGMFGSGGGGGGTGSTGPGYSFPHYGFPTYGGITFHPGTTKSNAGMKHGTMDTESKKDPEGCDKSDDKNTVNLFGKVIETTEQDQEPSEATVGNGEVTLTYATGTKEESAGVQDNLFLEKAMQLAKRHANALFDYAVTGDVKMLLAVQRHLTAVQDENGDSVLHLAIIHLHSQLVRDLLEVTSGLISDDIINMRNDLYQTPLHLAVITKQEDVVEDLLRAGADLSLLDRLGNSVLHLAAKEGHDKVLSILLKHKKAALLLDHPNGDGLNAIHLAMMSNSLPCLLLLVAAGADVNAQEQKSGRTALHLAVEHDNISLAGCLLLEGDAHVDSTTYDGTTPLHIAAGRGSTRLAALLKAAGADPLVENFEPLYDLDDSWENAGEDEGVVPGTTPLDMATSWQVFDILNGKPYEPEFTSDDLLAQGDMKQLAEDVKLQLYKLLEIPDPDKNWATLAQKLGLGILNNAFRLSPAPSKTLMDNYEVSGGTVRELVEALRQMGYTEAIEVIQAASSPVKTTSQAHSLPLSPASTRQQIDELRDSDSVCDSGVETSFRKLSFTESLTSGASLLTLNKMPHDYGQEGPLEGKIXXMAEDDPYLGRPEQMFHLDPSLTHTIFNPEVFQPQMALPTDGPYLQILEQPKQRGFRFRYVCEGPSHGGLPGASSEKNKKSYPQVKICNYVGPAKVIVQLVTNGKNIHLHAHSLVGKHCEDGICTVTAGPKDMVVGFANLGILHVTKKKVFETLEARMTEACIRGYNPGLLVHPDLAYLQAEGGGDRQLGDREKELIRQAALQQTKEMDLSVVRLMFTAFLPDSTGSFTRRLEPVVSDAIYDSKAPNASNLKIVRMDRTAGCVTGGEEIYLLCDKVQKDDIQIRFYEEEENGGVWEGFGDFSPTDVHRQFAIVFKTPKYKDINITKPASVFVQLRRKSDLETSEPKPFLYYPEIKDKEEVQRKRQKLMPNFSDSFGGGSGAGAGGGGMFGSGGGGGGTGSTGPGYSFPHYGFPTYGGITFHPGTTKSNAGMKHGTMDTESKKDPEGCDKSDDKNTVNLFGKVIETTEQDQEPSEATVGNGEVTLTYATGTKEESAGVQDNLFLEKAMQLAKRHANALFDYAVTGDVKMLLAVQRHLTAVQDENGDSVLHLAIIHLHSQLVRDLLEVTSGLISDDIINMRNDLYQTPLHLAVITKQEDVVEDLLRAGADLSLLDRLGNSVLHLAAKEGHDKVLSILLKHKKAALLLDHPNGDGLNAIHLAMMSNSLPCLLLLVAAGADVNAQEQKSGRTALHLAVEHDNISLAGCLLLEGDAHVDSTTYDGTTPLHIAAGRGSTRLAALLKAAGADPLVENFEPLYDLDDSWENAGEDEGVVPGTTPLDMATSWQVFDILNGKPYEPEFTSDDLLAQGDMKQLAEDVKLQLYKLLEIPDPDKNWATLAQKLGLGILNNAFRLSPAPSKTLMDNYEVSGGTVRELVEALRQMGYTEAIEVIQAASSPVKTTSQAHSLPLSPASTRQQIDELRDSDSVCDSGVETSFRKLSFTESLTSGASLLTLNKMPHDYGQEGPLEGKIXXMFHLDPSLTHTIFNPEVFQPQMALPTDGPYLQILEQPKQRGFRFRYVCEGPSHGGLPGASSEKNKKSYPQVKICNYVGPAKVIVQLVTNGKNIHLHAHSLVGKHCEDGICTVTAGPKDMVVGFANLGILHVTKKKVFETLEARMTEACIRGYNPGLLVHPDLAYLQAEGGGDRQLGDREKELIRQAALQQTKEMDLSVVRLMFTAFLPDSTGSFTRRLEPVVSDAIYDSKAPNASNLKIVRMDRTAGCVTGGEEIYLLCDKVQKDDIQIRFYEEEENGGVWEGFGDFSPTDVHRQFAIVFKTPKYKDINITKPASVFVQLRRKSDLETSEPKPFLYYPEIKDKEEVQRKRQKLMPNFSDSFGGGSGAGAGGGGMFGSGGGGGGTGSTGPGYSFPHYGFPTYGGITFHPGTTKSNAGMKHGTMDTESKKDPEGCDKSDDKNTVNLFGKVIETTEQDQEPSEATVGNGEVTLTYATGTKEESAGVQDNLFLEKAMQLAKRHANALFDYAVTGDVKMLLAVQRHLTAVQDENGDSVLHLAIIHLHSQLVRDLLEVTSGLISDDIINMRNDLYQTPLHLAVITKQEDVVEDLLRAGADLSLLDRLGNSVLHLAAKEGHDKVLSILLKHKKAALLLDHPNGDGLNAIHLAMMSNSLPCLLLLVAAGADVNAQEQKSGRTALHLAVEHDNISLAGCLLLEGDAHVDSTTYDGTTPLHIAAGRGSTRLAALLKAAGADPLVENFEPLYDLDDSWENAGEDEGVVPGTTPLDMATSWQVFDILNGKPYEPEFTSDDLLAQGDMKQLAEDVKLQLYKLLEIPDPDKNWATLAQKLGLGILNNAFRLSPAPSKTLMDNYEVSGGTVRELVEALRQMGYTEAIEVIQAASSPVKTTSQAHSLPLSPASTRQQIDELRDSDSVCDSGVETSFRKLSFTESLTSGASLLTLNKMPHDYGQEGPLEGKIXXMAEDDPYLGRPEQMFHLDPSLTHTIFNPEVFQPQMALPTADGPYLQILEQPKQRGFRFRYVCEGPSHGGLPGASSEKNKKSYPQVKICNYVGPAKVIVQLVTNGKNIHLHAHSLVGKHCEDGICTVTAGPKDMVVGFANLGILHVTKKKVFETLEARMTEACIRGYNPGLLVHPDLAYLQAEGGGDRQLGDREKELIRQAALQQTKEMDLSVVRLMFTAFLPDSTGSFTRRLEPVVSDAIYDSKAPNASNLKIVRMDRTAGCVTGGEEIYLLCDKVQKDDIQIRFYEEEENGGVWEGFGDFSPTDVHRQFAIVFKTPKYKDINITKPASVFVQLRRKSDLETSEPKPFLYYPEIKDKEEVQRKRQKLMPNFSDSFGGGSGAGAGGGGMFGSGGGGGGTGSTGPGYSFPHYGFPTYGGITFHPGTTKSNAGMKHGTMDTESKKDPEGCDKSDDKNTVNLFGKVIETTEQDQEPSEATVGNGEVTLTYATGTKEESAGVQDNLFLEKAMQLAKRHANALFDYAVTGDVKMLLAVQRHLTAVQDENGDSVLHLAIIHLHSQLVRDLLEVTSGLISDDIINMRNDLYQTPLHLAVITKQEDVVEDLLRAGADLSLLDRLGNSVLHLAAKEGHDKVLSILLKHKKAALLLDHPNGDGLNAIHLAMMSNSLPCLLLLVAAGADVNAQEQKSGRTALHLAVEHDNISLAGCLLLEGDAHVDSTTYDGTTPLHIAAGRGSTRLAALLKAAGADPLVENFEPLYDLDDSWENAGEDEGVVPGTTPLDMATSWQVFDILNGKPYEPEFTSDDLLAQGDMKQLAEDVKLQLYKLLEIPDPDKNWATLAQKLGLGILNNAFRLSPAPSKTLMDNYEVSGGTVRELVEALRQMGYTEAIEVIQAASSPVKTTSQAHSLPLSPASTRQQIDELRDSDSVCDSGVETSFRKLSFTESLTSGASLLTLNKMPHDYGQEGPLEGKIX"
## [6,] "MPRCPAGAMDEGPVDLRTRPKAAGLPGAALPLRKRPLRAPSPEPAAPRGAAGLVVPLDPLRGGCDLPAVPGPPHGLARPEALYYPGALLPLYPTRAMGSPFPLVNLPTPLYPMMCPMEHPLSADIAMATRADEDGDTPLHIAVVQGNLPAVHRLVNLFQQGGRELDIYNNLRQTPLHLAVITTLPSVVRLLVTAGASPMALDRHGQTAAHLACEHRSPTCLRALLDSAAPGTLDLEARNYDGLTALHVAVNTECQETVQLLLERGADIDAVDIKSGRSPLIHAVENNSLSMVQLLLQHGANVNAQMYSGSSALHSASGRGLLPLVRTLVRSGADSSLKNCHNDTPLMVARSRRVIDILRGKATRPASTSQPDPSPDRSANTSPESSSRLSSNGLLSASPSSSPSQSPPRDPPGFPMAPPNFFLPSPSPPAFLPFAGVLRGPGRPVPPSPAPGGSXXMAEDDPYLGRPEQMFHLDPSLTHTIFNPEVFQPQMALPTDGPYLQILEQPKQRGFRFRYVCEGPSHGGLPGASSEKNKKSYPQVKICNYVGPAKVIVQLVTNGKNIHLHAHSLVGKHCEDGICTVTAGPKDMVVGFANLGILHVTKKKVFETLEARMTEACIRGYNPGLLVHPDLAYLQAEGGGDRQLGDREKELIRQAALQQTKEMDLSVVRLMFTAFLPDSTGSFTRRLEPVVSDAIYDSKAPNASNLKIVRMDRTAGCVTGGEEIYLLCDKVQKDDIQIRFYEEEENGGVWEGFGDFSPTDVHRQFAIVFKTPKYKDINITKPASVFVQLRRKSDLETSEPKPFLYYPEIKDKEEVQRKRQKLMPNFSDSFGGGSGAGAGGGGMFGSGGGGGGTGSTGPGYSFPHYGFPTYGGITFHPGTTKSNAGMKHGTMDTESKKDPEGCDKSDDKNTVNLFGKVIETTEQDQEPSEATVGNGEVTLTYATGTKEESAGVQDNLFLEKAMQLAKRHANALFDYAVTGDVKMLLAVQRHLTAVQDENGDSVLHLAIIHLHSQLVRDLLEVTSGLISDDIINMRNDLYQTPLHLAVITKQEDVVEDLLRAGADLSLLDRLGNSVLHLAAKEGHDKVLSILLKHKKAALLLDHPNGDGLNAIHLAMMSNSLPCLLLLVAAGADVNAQEQKSGRTALHLAVEHDNISLAGCLLLEGDAHVDSTTYDGTTPLHIAAGRGSTRLAALLKAAGADPLVENFEPLYDLDDSWENAGEDEGVVPGTTPLDMATSWQVFDILNGKPYEPEFTSDDLLAQGDMKQLAEDVKLQLYKLLEIPDPDKNWATLAQKLGLGILNNAFRLSPAPSKTLMDNYEVSGGTVRELVEALRQMGYTEAIEVIQAASSPVKTTSQAHSLPLSPASTRQQIDELRDSDSVCDSGVETSFRKLSFTESLTSGASLLTLNKMPHDYGQEGPLEGKIXXMAEDDPYLGRPEQMFHLDPSLTHTIFNPEVFQPQMALPTADGPYLQILEQPKQRGFRFRYVCEGPSHGGLPGASSEKNKKSYPQVKICNYVGPAKVIVQLVTNGKNIHLHAHSLVGKHCEDGICTVTAGPKDMVVGFANLGILHVTKKKVFETLEARMTEACIRGYNPGLLVHPDLAYLQAEGGGDRQLGDREKELIRQAALQQTKEMDLSVVRLMFTAFLPDSTGSFTRRLEPVVSDAIYDSKAPNASNLKIVRMDRTAGCVTGGEEIYLLCDKVQKDDIQIRFYEEEENGGVWEGFGDFSPTDVHRQFAIVFKTPKYKDINITKPASVFVQLRRKSDLETSEPKPFLYYPEIKDKEEVQRKRQKLMPNFSDSFGGGSGAGAGGGGMFGSGGGGGGTGSTGPGYSFPHYGFPTYGGITFHPGTTKSNAGMKHGTMDTESKKDPEGCDKSDDKNTVNLFGKVIETTEQDQEPSEATVGNGEVTLTYATGTKEESAGVQDNLFLEKAMQLAKRHANALFDYAVTGDVKMLLAVQRHLTAVQDENGDSVLHLAIIHLHSQLVRDLLEVTSGLISDDIINMRNDLYQTPLHLAVITKQEDVVEDLLRAGADLSLLDRLGNSVLHLAAKEGHDKVLSILLKHKKAALLLDHPNGDGLNAIHLAMMSNSLPCLLLLVAAGADVNAQEQKSGRTALHLAVEHDNISLAGCLLLEGDAHVDSTTYDGTTPLHIAAGRGSTRLAALLKAAGADPLVENFEPLYDLDDSWENAGEDEGVVPGTTPLDMATSWQVFDILNGKPYEPEFTSDDLLAQGDMKQLAEDVKLQLYKLLEIPDPDKNWATLAQKLGLGILNNAFRLSPAPSKTLMDNYEVSGGTVRELVEALRQMGYTEAIEVIQAASSPVKTTSQAHSLPLSPASTRQQIDELRDSDSVCDSGVETSFRKLSFTESLTSGASLLTLNKMPHDYGQEGPLEGKIXXMAEDDPYLGRPEQMFHLDPSLTHTIFNPEVFQPQMALPTDGPYLQILEQPKQRGFRFRYVCEGPSHGGLPGASSEKNKKSYPQVKICNYVGPAKVIVQLVTNGKNIHLHAHSLVGKHCEDGICTVTAGPKDMVVGFANLGILHVTKKKVFETLEARMTEACIRGYNPGLLVHPDLAYLQAEGGGDRQLGDREKELIRQAALQQTKEMDLSVVRLMFTAFLPDSTGSFTRRLEPVVSDAIYDSKAPNASNLKIVRMDRTAGCVTGGEEIYLLCDKVQKDDIQIRFYEEEENGGVWEGFGDFSPTDVHRQFAIVFKTPKYKDINITKPASVFVQLRRKSDLETSEPKPFLYYPEIKDKEEVQRKRQKLMPNFSDSFGGGSGAGAGGGGMFGSGGGGGGTGSTGPGYSFPHYGFPTYGGITFHPGTTKSNAGMKHGTMDTESKKDPEGCDKSDDKNTVNLFGKVIETTEQDQEPSEATVGNGEVTLTYATGTKEESAGVQDNLFLEKAMQLAKRHANALFDYAVTGDVKMLLAVQRHLTAVQDENGDSVLHLAIIHLHSQLVRDLLEVTSGLISDDIINMRNDLYQTPLHLAVITKQEDVVEDLLRAGADLSLLDRLGNSVLHLAAKEGHDKVLSILLKHKKAALLLDHPNGDGLNAIHLAMMSNSLPCLLLLVAAGADVNAQEQKSGRTALHLAVEHDNISLAGCLLLEGDAHVDSTTYDGTTPLHIAAGRGSTRLAALLKAAGADPLVENFEPLYDLDDSWENAGEDEGVVPGTTPLDMATSWQVFDILNGKPYEPEFTSDDLLAQGDMKQLAEDVKLQLYKLLEIPDPDKNWATLAQKLGLGILNNAFRLSPAPSKTLMDNYEVSGGTVRELVEALRQMGYTEAIEVIQAASSPVKTTSQAHSLPLSPASTRQQIDELRDSDSVCDSGVETSFRKLSFTESLTSGASLLTLNKMPHDYGQEGPLEGKIXXMAEDDPYLGRPEQMFHLDPSLTHTIFNPEVFQPQMALPTADGPYLQILEQPKQRGFRFRYVCEGPSHGGLPGASSEKNKKSYPQVKICNYVGPAKVIVQLVTNGKNIHLHAHSLVGKHCEDGICTVTAGPKDMVVGFANLGILHVTKKKVFETLEARMTEACIRGYNPGLLVHPDLAYLQAEGGGDRQLGDREKELIRQAALQQTKEMDLSVVRLMFTAFLPDSTGSFTRRLEPVVSDAIYDSKAPNASNLKIVRMDRTAGCVTGGEEIYLLCDKVQKDDIQIRFYEEEENGGVWEGFGDFSPTDVHRQFAIVFKTPKYKDINITKPASVFVQLRRKSDLETSEPKPFLYYPEIKDKEEVQRKRQKLMPNFSDSFGGGSGAGAGGGGMFGSGGGGGGTGSTGPGYSFPHYGFPTYGGITFHPGTTKSNAGMKHGTMDTESKKDPEGCDKSDDKNTVNLFGKVIETTEQDQEPSEATVGNGEVTLTYATGTKEESAGVQDNLFLEKAMQLAKRHANALFDYAVTGDVKMLLAVQRHLTAVQDENGDSVLHLAIIHLHSQLVRDLLEVTSGLISDDIINMRNDLYQTPLHLAVITKQEDVVEDLLRAGADLSLLDRLGNSVLHLAAKEGHDKVLSILLKHKKAALLLDHPNGDGLNAIHLAMMSNSLPCLLLLVAAGADVNAQEQKSGRTALHLAVEHDNISLAGCLLLEGDAHVDSTTYDGTTPLHIAAGRGSTRLAALLKAAGADPLVENFEPLYDLDDSWENAGEDEGVVPGTTPLDMATSWQVFDILNGKPYEPEFTSDDLLAQGDMKQLAEDVKLQLYKLLEIPDPDKNWATLAQKLGLGILNNAFRLSPAPSKTLMDNYEVSGGTVRELVEALRQMGYTEAIEVIQAASSPVKTTSQAHSLPLSPASTRQQIDELRDSDSVCDSGVETSFRKLSFTESLTSGASLLTLNKMPHDYGQEGPLEGKIXXMAEDDPYLGRPEQMFHLDPSLTHTIFNPEVFQPQMALPTDGPYLQILEQPKQRGFRFRYVCEGPSHGGLPGASSEKNKKSYPQVKICNYVGPAKVIVQLVTNGKNIHLHAHSLVGKHCEDGICTVTAGPKDMVVGFANLGILHVTKKKVFETLEARMTEACIRGYNPGLLVHPDLAYLQAEGGGDRQLGDREKELIRQAALQQTKEMDLSVVRLMFTAFLPDSTGSFTRRLEPVVSDAIYDSKAPNASNLKIVRMDRTAGCVTGGEEIYLLCDKVQKDDIQIRFYEEEENGGVWEGFGDFSPTDVHRQFAIVFKTPKYKDINITKPASVFVQLRRKSDLETSEPKPFLYYPEIKDKEEVQRKRQKLMPNFSDSFGGGSGAGAGGGGMFGSGGGGGGTGSTGPGYSFPHYGFPTYGGITFHPGTTKSNAGMKHGTMDTESKKDPEGCDKSDDKNTVNLFGKVIETTEQDQEPSEATVGNGEVTLTYATGTKEESAGVQDNLFLEKAMQLAKRHANALFDYAVTGDVKMLLAVQRHLTAVQDENGDSVLHLAIIHLHSQLVRDLLEVTSGLISDDIINMRNDLYQTPLHLAVITKQEDVVEDLLRAGADLSLLDRLGNSVLHLAAKEGHDKVLSILLKHKKAALLLDHPNGDGLNAIHLAMMSNSLPCLLLLVAAGADVNAQEQKSGRTALHLAVEHDNISLAGCLLLEGDAHVDSTTYDGTTPLHIAAGRGSTRLAALLKAAGADPLVENFEPLYDLDDSWENAGEDEGVVPGTTPLDMATSWQVFDILNGKPYEPEFTSDDLLAQGDMKQLAEDVKLQLYKLLEIPDPDKNWATLAQKLGLGILNNAFRLSPAPSKTLMDNYEVSGGTVRELVEALRQMGYTEAIEVIQAASSPVKTTSQAHSLPLSPASTRQQIDELRDSDSVCDSGVETSFRKLSFTESLTSGASLLTLNKMPHDYGQEGPLEGKIXXMFHLDPSLTHTIFNPEVFQPQMALPTDGPYLQILEQPKQRGFRFRYVCEGPSHGGLPGASSEKNKKSYPQVKICNYVGPAKVIVQLVTNGKNIHLHAHSLVGKHCEDGICTVTAGPKDMVVGFANLGILHVTKKKVFETLEARMTEACIRGYNPGLLVHPDLAYLQAEGGGDRQLGDREKELIRQAALQQTKEMDLSVVRLMFTAFLPDSTGSFTRRLEPVVSDAIYDSKAPNASNLKIVRMDRTAGCVTGGEEIYLLCDKVQKDDIQIRFYEEEENGGVWEGFGDFSPTDVHRQFAIVFKTPKYKDINITKPASVFVQLRRKSDLETSEPKPFLYYPEIKDKEEVQRKRQKLMPNFSDSFGGGSGAGAGGGGMFGSGGGGGGTGSTGPGYSFPHYGFPTYGGITFHPGTTKSNAGMKHGTMDTESKKDPEGCDKSDDKNTVNLFGKVIETTEQDQEPSEATVGNGEVTLTYATGTKEESAGVQDNLFLEKAMQLAKRHANALFDYAVTGDVKMLLAVQRHLTAVQDENGDSVLHLAIIHLHSQLVRDLLEVTSGLISDDIINMRNDLYQTPLHLAVITKQEDVVEDLLRAGADLSLLDRLGNSVLHLAAKEGHDKVLSILLKHKKAALLLDHPNGDGLNAIHLAMMSNSLPCLLLLVAAGADVNAQEQKSGRTALHLAVEHDNISLAGCLLLEGDAHVDSTTYDGTTPLHIAAGRGSTRLAALLKAAGADPLVENFEPLYDLDDSWENAGEDEGVVPGTTPLDMATSWQVFDILNGKPYEPEFTSDDLLAQGDMKQLAEDVKLQLYKLLEIPDPDKNWATLAQKLGLGILNNAFRLSPAPSKTLMDNYEVSGGTVRELVEALRQMGYTEAIEVIQAASSPVKTTSQAHSLPLSPASTRQQIDELRDSDSVCDSGVETSFRKLSFTESLTSGASLLTLNKMPHDYGQEGPLEGKIXXMAEDDPYLGRPEQMFHLDPSLTHTIFNPEVFQPQMALPTADGPYLQILEQPKQRGFRFRYVCEGPSHGGLPGASSEKNKKSYPQVKICNYVGPAKVIVQLVTNGKNIHLHAHSLVGKHCEDGICTVTAGPKDMVVGFANLGILHVTKKKVFETLEARMTEACIRGYNPGLLVHPDLAYLQAEGGGDRQLGDREKELIRQAALQQTKEMDLSVVRLMFTAFLPDSTGSFTRRLEPVVSDAIYDSKAPNASNLKIVRMDRTAGCVTGGEEIYLLCDKVQKDDIQIRFYEEEENGGVWEGFGDFSPTDVHRQFAIVFKTPKYKDINITKPASVFVQLRRKSDLETSEPKPFLYYPEIKDKEEVQRKRQKLMPNFSDSFGGGSGAGAGGGGMFGSGGGGGGTGSTGPGYSFPHYGFPTYGGITFHPGTTKSNAGMKHGTMDTESKKDPEGCDKSDDKNTVNLFGKVIETTEQDQEPSEATVGNGEVTLTYATGTKEESAGVQDNLFLEKAMQLAKRHANALFDYAVTGDVKMLLAVQRHLTAVQDENGDSVLHLAIIHLHSQLVRDLLEVTSGLISDDIINMRNDLYQTPLHLAVITKQEDVVEDLLRAGADLSLLDRLGNSVLHLAAKEGHDKVLSILLKHKKAALLLDHPNGDGLNAIHLAMMSNSLPCLLLLVAAGADVNAQEQKSGRTALHLAVEHDNISLAGCLLLEGDAHVDSTTYDGTTPLHIAAGRGSTRLAALLKAAGADPLVENFEPLYDLDDSWENAGEDEGVVPGTTPLDMATSWQVFDILNGKPYEPEFTSDDLLAQGDMKQLAEDVKLQLYKLLEIPDPDKNWATLAQKLGLGILNNAFRLSPAPSKTLMDNYEVSGGTVRELVEALRQMGYTEAIEVIQAASSPVKTTSQAHSLPLSPASTRQQIDELRDSDSVCDSGVETSFRKLSFTESLTSGASLLTLNKMPHDYGQEGPLEGKIX"
##      Mutant_Peptide                Total_RNA Tumor_RNA_Ratio Tumor_RNA
## [1,] "MAEDDPYLGRPEKMFHLDPSLTHTIFN" "NA"      "NA"            "NA"     
## [2,] "MAEDDPYLGRPEKMFHLDPSLTHTIFN" "NA"      "NA"            "NA"     
## [3,] "MAEDDPYLGRPEKMFHLDPSLTHTIFN" "NA"      "NA"            "NA"     
## [4,] "MAEDDPYLGRPEKMFHLDPSLTHTIFN" "NA"      "NA"            "NA"     
## [5,] "MAEDDPYLGRPEKMFHLDPSLTHTIFN" "NA"      "NA"            "NA"     
## [6,] "MAEDDPYLGRPEKMFHLDPSLTHTIFN" "NA"      "NA"            "NA"     
##      Tumor_RNA_based_on_DNA MutRatio MutRatio_Min MutRatio_Max P_I P_R P 
## [1,] "NA"                   "NA"     "NA"         "NA"         NA  NA  NA
## [2,] "NA"                   "NA"     "NA"         "NA"         NA  NA  NA
## [3,] "NA"                   "NA"     "NA"         "NA"         NA  NA  NA
## [4,] "NA"                   "NA"     "NA"         "NA"         NA  NA  NA
## [5,] "NA"                   "NA"     "NA"         "NA"         NA  NA  NA
## [6,] "NA"                   "NA"     "NA"         "NA"         NA  NA  NA
```


``` r
  print(Export_Summary_Fragments(Result_HLA2_Seq, Mut_Rank_th = 5))
```

```
## [[1]]
##                                                     atggcagaag
## Num_Peptide_Per_Grp                                     13.000
## Num_Cond_Peptide_Per_Grp                                13.000
## Num_Rest_Peptide_Per_Grp                                 3.000
## Num_Rest_Peptide_Per_Grp / Num_Cond_Peptide_Per_Grp      0.231
## -logP                                                    0.621
## 
## [[2]]
##                                                         
## Num_Peptide_Per_NM                                13.000
## Num_Cond_Peptide_Per_NM                           13.000
## Num_Rest_Peptide_Per_NM                            3.000
## Num_Rest_Peptide_Per_NM / Num_Cond_Peptide_Per_NM  0.231
## -logP                                              0.621
## 
## [[3]]
##                                                     MAEDDPYLGRPEKMFHLDPSLTHTIFN-0_atggcagaag-0_1
## Num_Peptide_Per_Pep                                                                       13.000
## Num_Cond_Peptide_Per_Pep                                                                  13.000
## Num_Rest_Peptide_Per_Pep                                                                   3.000
## Num_Rest_Peptide_Per_Pep / Num_Cond_Peptide_Per_Pep                                        0.231
## -logP                                                                                      0.621
```
