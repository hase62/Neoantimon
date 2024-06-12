#'Execute Sample Analysis
#'
#'@return void
#'
#'@export
TestAnalysis<-function(){
  #library(devtools);
  #install_github('hase62/Neoantimon');
  library(Neoantimon);
  #library(biomaRt)

  data("sample_vcf.annovar")
  data("sample_vcf.vep")
  data("sample_hla_table_c1")
  data("sample_refFlat.grch37")
  data("sample_refMrna.grch37.fa")
  data("sample_result_SNV_CLASS1_ALL")

  for(r_code in dir("./../R/", ".R")) source(paste("./../R/", r_code, sep = ""))

  MainSNVClass1(input_annovar_format_file = sample_vcf.annovar,
                hla_types = sample_hla_table_c1[1,-1],
                refflat_file = sample_refFlat.grch37,
                refmrna_file = sample_refMrna.grch37.fa,
                netMHCpan_dir = NA)

  write.table(file = "result.ID.SNV1/data.ID_SNV.peptide.SNV_CLASS1_ALL.txt",
              x = sample_result_SNV_CLASS1_ALL[grep("0_DHX15", sample_result_SNV_CLASS1_ALL$Gene), ],
              row.names = FALSE, quote = FALSE, sep = "\t")

  data("sample_vcf.vep")
  print(sample_vcf.vep, row.names = FALSE)

  data("sample_vcf.snps")
  print(sample_vcf.snps, row.row.names = FALSE)

  Result_HLA1_SNV <- MainSNVClass1(input_annovar_format_file = "data/sample_vcf.annovar.txt",
                                   file_name_in_hla_table = "sample",
                                   hla_file = "data/sample_hla_table_c1.txt",
                                   refflat_file  = "refFlat.grch37.txt",
                                   refmrna_file = "refMrna.grch37.fa",
                                   rnaexp_file = "data/sample_rna_exp.txt",
                                   netMHCpan_dir = "netMHCpan-4.1/netMHCpan",
                                   depth_tumor_column = 12,
                                   depth_normal_column = 14,
                                   multiple_variants = TRUE,
                                   SNPs = "data/sample_vcf.snps.txt")
  Result_HLA1_SNV_1 <- CalculatePriorityScores(result = Result_HLA1_SNV[[1]], useRNAvaf = FALSE)
  Result_HLA1_SNV_2 <- CalculatePriorityScores(result = Result_HLA1_SNV[[2]], useRNAvaf = FALSE)
  Export_Summary_SNV(Input = Result_HLA1_SNV_1[[1]], Mut_Rank_th = 0.05, Wt_Rank_th = 0.05)
  Export_Summary_SNV(Input = Result_HLA1_SNV_1[[2]], Mut_Rank_th = 0.05, Wt_Rank_th = 0.05)
  write.table(Exp, "Exp.txt", row.names = TRUE, col.names = FALSE, quote = FALSE, sep = "\t")

  Result_HLA1_SNV_vep <- MainSNVClass1(input_vep_format_file = "data/sample_vcf.vep.txt",
                                   file_name_in_hla_table = "sample",
                                   hla_file = "data/sample_hla_table_c1.txt",
                                   refflat_file  = "refFlat.grch37.txt",
                                   refmrna_file = "refMrna.grch37.fa",
                                   rnaexp_file = "data/sample_rna_exp.txt",
                                   netMHCpan_dir = "netMHCpan-4.1/netMHCpan",
                                   multiple_variants = FALSE)
  Result_HLA1_vep_SNV <- CalculatePriorityScores(result = Result_HLA1_SNV_vep, useRNAvaf = FALSE)
  Export_Summary_SNV(Input = Result_HLA1_vep_SNV, Mut_Rank_th = 0.05, Wt_Rank_th = 0.05)

  Result_HLA2_SNV <- MainSNVClass2(input_annovar_format_file = "data/sample_vcf.annovar.txt",
                                   file_name_in_hla_table = "sample",
                                   hla_file = "data/sample_hla_table_c2.txt",
                                   refflat_file  = "refFlat.grch37.txt",
                                   refmrna_file = "refMrna.grch37.fa",
                                   rnaexp_file = "data/sample_rna_exp.txt",
                                   netMHCIIpan_dir = "netMHCIIpan-4.2/netMHCIIpan",
                                   depth_tumor_column = 12,
                                   depth_normal_column = 14,
                                   SNPs = "data/sample_vcf.snps.vcf",
                                   multiple_variants = TRUE,
                                   MHCflurry = "~/opt/anaconda3/bin/mhctools")
  Result_HLA2_SNV <- CalculatePriorityScores(result = Result_HLA2_SNV, useRNAvaf = FALSE)

  Result_HLA1_INDEL <- MainINDELClass1(input_annovar_format_file = "data/sample_vcf.annovar.txt",
                                       file_name_in_hla_table = "sample",
                                       hla_file = "data/sample_hla_table_c1.txt",
                                       refflat_file  = "refFlat.grch37.txt",
                                       refmrna_file = "refMrna.grch37.fa",
                                       rnaexp_file = "data/sample_rna_exp.txt",
                                       netMHCpan_dir = "netMHCpan-4.1/netMHCpan",
                                       depth_tumor_column = 12,
                                       depth_normal_column = 14,
                                       multiple_variants = TRUE,
                                       base_0 = FALSE)
  Result_HLA1_INDEL_1 <- CalculatePriorityScores(result = Result_HLA1_INDEL[[1]], useRNAvaf = FALSE)
  Result_HLA1_INDEL_2 <- CalculatePriorityScores(result = Result_HLA1_INDEL[[2]], useRNAvaf = FALSE)
  Export_Summary_IndelSV(Input = Result_HLA1_INDEL, Mut_Rank_th = 0.05)

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
  Result_HLA2_INDEL <- CalculatePriorityScores(result = Result_HLA2_INDEL, useRNAvaf = FALSE)

  Result_HLA1_SV <- MainSVFUSIONClass1(input_file = "data/sample_sv_bnd.txt",
                                       file_name_in_hla_table = "sample",
                                       hla_file = "data/sample_hla_table_c1.txt",
                                       refflat_file  = "refFlat.grch37.txt",
                                       refmrna_file = "refMrna.grch37.fa",
                                       rnaexp_file = "data/sample_rna_exp.txt",
                                       netMHCpan_dir = "netMHCpan-4.0/netMHCpan",
                                       mutation_alt_bnd_column = 5,
                                       gene_symbol_column = 7,
                                       mate_id_column = 8)
  Result_HLA1_SV <- CalculatePriorityScores(result = Result_HLA1_SV, useRNAvaf = FALSE)

  Result_HLA2_SV <- MainSVFUSIONClass2(input_file = "data/sample_sv_bnd.txt",
                                       file_name_in_hla_table = "sample",
                                       hla_file = "data/sample_hla_table_c2.txt",
                                       refflat_file  = "refFlat.grch37.txt",
                                       refmrna_file = "refMrna.grch37.fa",
                                       rnaexp_file = "data/sample_rna_exp.txt",
                                       netMHCIIpan_dir = "netMHCIIpan-3.2/netMHCIIpan",
                                       mutation_alt_bnd_column = 5,
                                       gene_symbol_column = 7,
                                       mate_id_column = 8)
  Result_HLA2_SV <- CalculatePriorityScores(result = Result_HLA2_SV, useRNAvaf = FALSE)

  Result_HLA1_Seq <- MainSeqFragmentClass1(input_sequence = "atggcagaagatgatccatatttgggaaggcctgaaaaaatgtttcatttggatccttctttgactcatacaatatttaatc",
                                           file_name_in_hla_table = "sample",
                                           hla_file = "data/sample_hla_table_c1.txt",
                                           hmdir = getwd(),
                                           job_id = "NO_job_id",
                                           refflat_file  = "refFlat.grch37.txt",
                                           refmrna_file = "refMrna.grch37.fa",
                                           netMHCpan_dir = "netMHCpan-4.0/netMHCpan",
                                           reference_nm_id = c("NM_003998", "NM_001165412"))
  Result_HLA1_Seq <- CalculatePriorityScores(result = Result_HLA1_Seq, useRNAvaf = FALSE)

  Result_HLA2_Seq <- MainSeqFragmentClass2(input_sequence = "atggcgagcagcgagcagccagcagcaatgcaggcagcagaagatgatccatatttgggaaggcctgaacaaatgtttcatttgatccttctttgactcatacaatatttaatc",
                                           file_name_in_hla_table = "sample",
                                           hla_file = "data/sample_hla_table_c2.txt",
                                           hmdir = getwd(),
                                           job_id = "NO_job_id",
                                           refflat_file  = "refFlat.grch37.txt",
                                           refmrna_file = "refMrna.grch37.fa",
                                           netMHCIIpan_dir = "netMHCIIpan-3.2/netMHCIIpan",
                                           reference_gene_symbol = c("NFKB1", "BCL3"))
  Result_HLA2_Seq <- CalculatePriorityScores(result = Result_HLA2_Seq, useRNAvaf = FALSE)
}
