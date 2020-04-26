#'Execute Sample Analysis
#'
#'@return void
#'
#'@export
TestAnalysis<-function(){
  #library(devtools);
  #install_github('hase62/Neoantimon');
  #library(Neoantimon);
  #library(biomaRt)

  data("sample_vcf.annovar")
  data("sample_vcf.vep")
  data("sample_hla_table_c1")
  data("sample_refFlat.grch37")
  data("sample_refMrna.grch37.fa")
  data("sample_result_SNV_CLASS1_ALL")

  MainSNVClass1(input_annovar_format_file = sample_vcf.annovar,
                hla_types = sample_hla_table_c1[1,-1],
                refflat_file = sample_refFlat.grch37,
                refmrna_file = sample_refMrna.grch37.fa,
                netMHCpan_dir = NA)

  MainSNVClass1(input_vep_format_file = sample_vcf.vep,
                hla_types = sample_hla_table_c1[1,-1],
                refflat_file = sample_refFlat.grch37,
                refmrna_file = sample_refMrna.grch37.fa,
                netMHCpan_dir = NA)

  write.table(file = "result.ID.SNV1/data.ID_SNV.peptide.SNV_CLASS1_ALL.txt",
              x = sample_result_SNV_CLASS1_ALL[grep("0_DHX15", sample_result_SNV_CLASS1_ALL$Gene), ],
              row.names = FALSE, quote = FALSE, sep = "\t")

  data("sample_vcf.vep")
  print(sample_vcf.vep, row.names = FALSE)

  data("sample.snps.vcf")
  print(sample.snps.vcf, row.row.names = FALSE)

  Result_HLA1_INDEL <- MainSNVClass1(input_vep_format_file = "data/sample_vcf.vep.txt",
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

  Result_HLA1_Seq <- MainSeqFragmentClass1(input_sequence = "atggcagaagatgatccatatttgggaaggcctgaaaaaatgtttcatttggatccttctttgactcatacaatatttaatc",
                                           file_name_in_hla_table = "sample",
                                           hla_file = "data/sample_hla_table_c1.txt",
                                           hmdir = getwd(),
                                           job_id = "NO_job_id",
                                           refflat_file  = "refFlat.grch37.txt",
                                           refmrna_file = "refMrna.grch37.fa",
                                           netMHCpan_dir = "netMHCpan-4.0/netMHCpan",
                                           reference_nm_id = c("NM_003998", "NM_001165412"))

  Result_HLA2_Seq <- MainSeqFragmentClass2(input_sequence = "atggcagaagatgatccatatttgggaaggcctgaacaaatgtttcatttgatccttctttgactcatacaatatttaatc",
                                           file_name_in_hla_table = "sample",
                                           hla_file = "data/sample_hla_table_c2.txt",
                                           hmdir = getwd(),
                                           job_id = "NO_job_id",
                                           refflat_file  = "refFlat.grch37.txt",
                                           refmrna_file = "refMrna.grch37.fa",
                                           netMHCIIpan_dir = "netMHCIIpan-3.2/netMHCIIpan",
                                           reference_gene_symbol = c("NFKB1", "BCL3"))
}
