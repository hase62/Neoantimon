#'Execute Sample Analysis
#'
#'@return void
#'
#'@export
TestAnalysis<-function(){
  install.packages('devtools');
  library(devtools);
  install_github('hase62/Neoantimon');
  library(Neoantimon);
  #ls(all=T)
  #character(0)
  #rm(list=ls(all=TRUE))

  #for(rc in list.files("/Users/takaorihasegawa/Git/Neoantimon/R")) source(paste("/Users/takaorihasegawa/Git/Neoantimon/R", rc, sep = "/"))

  print("Please Install NetMHCpan and NetMHCIIpan if you did not do it.")
  print("Please Download refFlat and refmRNA files if you did not do it.")

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
  print(head(Result_HLA1_SNV))
  print(Export_Summary_SNV(Input = Result_HLA1_SNV, Mut_IC50_th = 500, Wt_IC50_th = 500))

  Result_HLA1_SNV <- MainSNVClass1(input_annovar_format_file =  "data/sample.vep.txt.annovar_format.txt",
                                   file_name_in_hla_table = "sample",
                                   hla_file = "data/sample_hla_table_c1.txt",
                                   refflat_file  = "lib/refFlat.grch38.txt",
                                   refmrna_file = "lib/refMrna.grch38.fa",
                                   rnaexp_file = "data/sample_rna_exp.txt",
                                   netMHCpan_dir = "lib/netMHCpan-4.0/netMHCpan",
                                   SNPs = "lib/sample.snps.vcf",
                                   multiple_variants = TRUE,
                                   MHCflurry = "mhctools")
  print(head(Result_HLA1_SNV))
  print(Export_Summary_SNV(Input = Result_HLA1_SNV, Mut_IC50_th = 500, Wt_IC50_th = 500))

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
  print(head(Result_HLA2_SNV))
  print(Export_Summary_SNV(Input = Result_HLA2_SNV, Mut_IC50_th = 500, Wt_IC50_th = 500))

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
  print(head(Result_HLA1_INDEL))
  print(Export_Summary_IndelSV(Input = Result_HLA1_INDEL, Mut_IC50_th = 500))
  print(Export_Summary_IndelSV_perFragments(Input = Result_HLA1_INDEL, Mut_IC50_th = 500))

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
  print(head(Result_HLA2_INDEL))

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
  print(head(Result_HLA1_SV))
  print(Export_Summary_IndelSV(Result_HLA1_SV, Mut_IC50_th = 500))
  print(Export_Summary_IndelSV_perFragments(Result_HLA1_SV, Mut_IC50_th = 500))

  Result_HLA2_SV <- MainSVFUSIONClass2(input_file = "data/sample_sv_bnd.txt",
                                       file_name_in_hla_table = "sample",
                                       hla_file = "data/sample_hla_table_c2.txt",
                                       refflat_file  = "lib/refFlat.txt",
                                       refmrna_file = "lib/refMrna.fa",
                                       rnaexp_file = "data/sample_rna_exp.txt",
                                       netMHCIIpan_dir = "netMHCIIpan-3.2/netMHCIIpan",
                                       refdna_file = "lib/GRCh37.fa",
                                       mutation_alt_bnd_column = 5,
                                       gene_symbol_column = 7,
                                       mate_id_column = 8)
  print(head(Result_HLA2_SV))

  Result_HLA1_Seq <- MainSeqFragmentClass1(input_sequence = "atggcagaagatgatccatatttgggaaggcctgaaaaaatgtttcatttggatccttctttgactcatacaatatttaatc",
                                           file_name_in_hla_table = "sample",
                                           hla_file = "data/sample_hla_table_c1.txt",
                                           hmdir = getwd(),
                                           job_id = "NO_job_id",
                                           refflat_file  = "lib/refFlat.txt",
                                           refmrna_file = "lib/refMrna.fa",
                                           netMHCpan_dir = "lib/netMHCpan-4.0/netMHCpan",
                                           reference_nm_id = c("NM_003998", "NM_001165412"))
  print(head(Result_HLA1_Seq))
  print(Export_Summary_Fragments(Result_HLA1_Seq, Mut_IC50_th = 500))

  Result_HLA2_Seq <- MainSeqFragmentClass2(input_sequence = "atggcagaagatgatccatatttgggaaggcctgaacaaatgtttcatttgatccttctttgactcatacaatatttaatc",
                                           file_name_in_hla_table = "sample",
                                           hla_file = "data/sample_hla_table_c2.txt",
                                           hmdir = getwd(),
                                           job_id = "NO_job_id",
                                           refflat_file  = "lib/refFlat.txt",
                                           refmrna_file = "lib/refMrna.fa",
                                           netMHCIIpan_dir = "lib/netMHCIIpan-3.1/netMHCIIpan",
                                           reference_gene_symbol = c("NFKB1", "BCL3")
  )
  print(head(Result_HLA2_Seq))
}
