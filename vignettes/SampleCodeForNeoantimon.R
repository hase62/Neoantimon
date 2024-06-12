## ----Preparation--------------------------------------------------------------
#install.packages('devtools');
library(devtools);
install_github('hase62/Neoantimon');
library(Neoantimon);
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("biomaRt")
library(biomaRt)

## ----Get SNV Sample 1 in Test Analysis----------------------------------------
Result_HLA1_SNV <- MainSNVClass1(input_annovar_format_file = "data/sample_vcf.annovar.txt",
                                   file_name_in_hla_table = "sample",
                                   hla_file = "data/sample_hla_table_c1.txt",
                                   refflat_file  = "refFlat.grch37.txt",
                                   refmrna_file = "refMrna.grch37.fa",
                                   rnaexp_file = "data/sample_rna_exp.txt",
                                   netMHCpan_dir = "netMHCpan-4.1/netMHCpan",
                                   depth_tumor_column = 12,
                                   depth_normal_column = 14,
                                   SNPs = "data/sample_vcf.snps.vcf",
                                   multiple_variants = TRUE)
  Result_HLA1_SNV_1 <- CalculatePriorityScores(result = Result_HLA1_SNV[[1]], useRNAvaf = FALSE)
  Result_HLA1_SNV_2 <- CalculatePriorityScores(result = Result_HLA1_SNV[[2]], useRNAvaf = FALSE)
  print(head(Result_HLA1_SNV_1))

## ----Get SNV Summary 1 in Test Analysis---------------------------------------
  print(Export_Summary_SNV(Input = Result_HLA1_SNV_1, Mut_Rank_th = 5, Wt_Rank_th = 5))

