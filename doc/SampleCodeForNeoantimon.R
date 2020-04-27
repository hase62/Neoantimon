## ----Preparation--------------------------------------------------------------
#install.packages('devtools');
library(devtools);
install_github('hase62/Neoantimon');
library(Neoantimon);

## ----Get SNV Sample 1 in Test Analysis----------------------------------------
Result_HLA1_SNV <- MainSNVClass1(input_annovar_format_file = "data/sample_vcf.annovar.txt",
                                   file_name_in_hla_table = "sample",
                                   hla_file = "data/sample_hla_table_c1.txt",
                                   refflat_file  = "refFlat.grch37.txt",
                                   refmrna_file = "refMrna.grch37.fa",
                                   rnaexp_file = "data/sample_rna_exp.txt",
                                   netMHCpan_dir = "netMHCpan-4.0/netMHCpan",
                                   depth_tumor_column = 12,
                                   depth_normal_column = 14,
                                   SNPs = "data/sample_vcf.snps.vcf",
                                   multiple_variants = TRUE,
                                   MHCflurry = "~/opt/anaconda3/bin/mhctools")
  print(head(Result_HLA1_SNV[[1]]))

## ----Get SNV Summary 1 in Test Analysis---------------------------------------
  print(Export_Summary_SNV(Input = Result_HLA1_SNV[[1]], Mut_IC50_th = 500, Wt_IC50_th = 500))

## ----Get SNV Sample 2 in Test Analysis----------------------------------------
 Result_HLA2_SNV <- MainSNVClass2(input_annovar_format_file = "data/sample_vcf.annovar.txt",
                                   file_name_in_hla_table = "sample",
                                   hla_file = "data/sample_hla_table_c2.txt",
                                   refflat_file  = "refFlat.grch37.txt",
                                   refmrna_file = "refMrna.grch37.fa",
                                   rnaexp_file = "data/sample_rna_exp.txt",
                                   netMHCIIpan_dir = "netMHCIIpan-3.2/netMHCIIpan",
                                   depth_tumor_column = 12,
                                   depth_normal_column = 14,
                                   SNPs = "data/sample_vcf.snps.vcf",
                                   multiple_variants = TRUE)
  print(head(Result_HLA2_SNV))

## ----Get SNV Summary 2 in Test Analysis---------------------------------------
  print(Export_Summary_SNV(Input = Result_HLA2_SNV, Mut_IC50_th = 500, Wt_IC50_th = 500))

## ----Get INDEL Sample 1 in Test Analysis--------------------------------------
Result_HLA1_INDEL <- MainINDELClass1(input_annovar_format_file = "data/sample_vcf.annovar.txt",
                                       file_name_in_hla_table = "sample",
                                       hla_file = "data/sample_hla_table_c1.txt",
                                       refflat_file  = "refFlat.grch37.txt",
                                       refmrna_file = "refMrna.grch37.fa",
                                       rnaexp_file = "data/sample_rna_exp.txt",
                                       netMHCpan_dir = "netMHCpan-4.0/netMHCpan",
                                       depth_tumor_column = 12,
                                       depth_normal_column = 14,
                                       SNPs = "data/sample_vcf.snps.vcf",
                                       multiple_variants = TRUE,
                                       MHCflurry = "~/opt/anaconda3/bin/mhctools")
  print(head(Result_HLA1_INDEL[[1]]))

## ----Get INDEL Summary 1-1 in Test Analysis-----------------------------------
  print(Export_Summary_IndelSV(Input = Result_HLA1_INDEL[[1]], Mut_IC50_th = 500))

## ----Get INDEL Summary 1-2 in Test Analysis-----------------------------------
  print(Export_Summary_IndelSV_perFragments(Input = Result_HLA1_INDEL[[1]], Mut_IC50_th = 500))

## ----Get INDEL Sample 2 in Test Analysis--------------------------------------
  Result_HLA2_INDEL <- MainINDELClass2(input_annovar_format_file = "data/sample_vcf.annovar.txt",
                                       file_name_in_hla_table = "sample",
                                       hla_file = "data/sample_hla_table_c2.txt",
                                       refflat_file  = "refFlat.grch37.txt",
                                       refmrna_file = "refMrna.grch37.fa",
                                       rnaexp_file = "data/sample_rna_exp.txt",
                                       netMHCIIpan_dir = "netMHCIIpan-3.2/netMHCIIpan",
                                       depth_tumor_column = 12,
                                       depth_normal_column = 14,
                                       SNPs = "data/sample_vcf.snps.vcf",
                                       multiple_variants = TRUE)
  print(head(Result_HLA2_INDEL))

## ----Get INDEL Summary 2-1 in Test Analysis-----------------------------------
  print(Export_Summary_IndelSV(Input = Result_HLA2_INDEL, Mut_IC50_th = 500))

## ----Get INDEL Summary 2-2 in Test Analysis-----------------------------------
  print(Export_Summary_IndelSV_perFragments(Input = Result_HLA2_INDEL, Mut_IC50_th = 500))

## ----Get SV Sample 1 in Test Analysis-----------------------------------------
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
  print(head(Result_HLA1_SV))

## ----Get SV Summary 1-1 in Test Analysis--------------------------------------
  print(Export_Summary_IndelSV(Result_HLA1_SV, Mut_IC50_th = 500))

## ----Get SV Summary 1-2 in Test Analysis--------------------------------------
  print(Export_Summary_IndelSV_perFragments(Result_HLA1_SV, Mut_IC50_th = 500))

## ----Get SV Sample 2 in Test Analysis-----------------------------------------
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
  print(head(Result_HLA2_SV))

## ----Get SV Summary 2-1 in Test Analysis--------------------------------------
  print(Export_Summary_IndelSV(Result_HLA2_SV, Mut_IC50_th = 500))

## ----Get SV Summary 2-2 in Test Analysis--------------------------------------
  print(Export_Summary_IndelSV_perFragments(Result_HLA2_SV, Mut_IC50_th = 500))

## ----Get Fragment Sample 1 in Test Analysis-----------------------------------
   Result_HLA1_Seq <- MainSeqFragmentClass1(input_sequence = "atggcagaagatgatccatatttgggaaggcctgaaaaaatgtttcatttggatccttctttgactcatacaatatttaatc",
                                           file_name_in_hla_table = "sample",
                                           hla_file = "data/sample_hla_table_c1.txt",
                                           hmdir = getwd(),
                                           job_id = "NO_job_id",
                                           refflat_file  = "refFlat.grch37.txt",
                                           refmrna_file = "refMrna.grch37.fa",
                                           netMHCpan_dir = "netMHCpan-4.0/netMHCpan",
                                           reference_nm_id = c("NM_003998", "NM_001165412"))
print(head(Result_HLA1_Seq))

## ----Get Fragment Summary 1 in Test Analysis----------------------------------
  print(Export_Summary_Fragments(Result_HLA1_Seq, Mut_IC50_th = 500))

## ----Get Fragment Sample 2 in Test Analysis-----------------------------------
  Result_HLA2_Seq <- MainSeqFragmentClass2(input_sequence = "atggcagaagatgatccatatttgggaaggcctgaacaaatgtttcatttgatccttctttgactcatacaatatttaatc",
                                           file_name_in_hla_table = "sample",
                                           hla_file = "data/sample_hla_table_c2.txt",
                                           hmdir = getwd(),
                                           job_id = "NO_job_id",
                                           refflat_file  = "refFlat.grch37.txt",
                                           refmrna_file = "refMrna.grch37.fa",
                                           netMHCIIpan_dir = "netMHCIIpan-3.2/netMHCIIpan",
                                           reference_gene_symbol = c("NFKB1", "BCL3"))
print(head(Result_HLA2_Seq))

## ----Get Fragment Summary 2 in Test Analysis----------------------------------
  print(Export_Summary_Fragments(Result_HLA2_Seq, Mut_IC50_th = 500))

