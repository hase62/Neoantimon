#'Execute Sample Analysis
#'
#'@return void
#'
#'@export
TestAnalysis<-function(){
  print("Please Install NetMHCpan and NetMHCIIpan if you did not do it.")
  InstallSampleFiles(export_dir = "lib")
  InstallSamtools(url = "https://github.com/hase62/Neoantimon/raw/master/lib/samtools-0.1.19.tar.bz2",
                  export_dir = "lib")
  InstallRefFlat(url = "http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/refFlat.txt.gz",
                 export_dir = "lib")
  InstallRefMrnaFile(url = "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/refMrna.fa.gz",
                     export_dir = "lib")

  Result_HLA1_SNV <- MainSNVClass1(input_file = "lib/data/sample_vcf.txt",
                                   file_name_in_hla_table = "sample",
                                   hla_file = "lib/data/sample_hla_table_c1.txt",
                                   refflat_file  = "lib/refFlat.txt",
                                   refmrna_file = "lib/refMrna.fa",
                                   rnaexp_file = "lib/data/sample_rna_exp.txt",
                                   netMHCpan_dir = "lib/netMHCpan-4.0/netMHCpan",
                                   nm_id_column = 10,
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
                                   nm_id_column = 10,
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
                                       nm_id_column = 10,
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
                                       nm_id_column = 10,
                                       depth_tumor_column = 12,
                                       depth_normal_column = 14
  )

  Result_HLA1_SV <- MainSVFUSIONClass1(input_file = "lib/data/sample_sv_bnd.txt",
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

  Result_HLA1_Seq <- MainSeqFragmentClass1(input_sequence = "atggcagaagatgatccatatttgggaaggcctgaaaaaatgtttcatttggatccttctttgactcatacaatatttaatc",
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

  Result_HLA2_Seq <- MainSeqFragmentClass2(input_sequence = "aaatggcagaagatgatccatatttgggaaggcctgaacaaatgtttcatttgatccttctttgactcatacaatatttaatc",
                                           file_name_in_hla_table = "sample",
                                           hla_file = "lib/data/sample_hla_table_c2.txt",
                                           hmdir = getwd(),
                                           job_id = "NO_job_id",
                                           refflat_file  = "lib/refFlat.txt",
                                           refmrna_file = "lib/refMrna.fa",
                                           netMHCIIpan_dir = "lib/netMHCIIpan-3.1/netMHCIIpan",
                                           gene_symbol = c("NFKB1", "BCL3"),
                                           reading_frame = 2
  )
}
