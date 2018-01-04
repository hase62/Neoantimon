#'Execute Sample Analysis
#'
#'@return void
#'
#'@export
TestAnalysis<-function(){
  print("Install NetMHCpan and NetMHCIIpan.")
  InstallSampleFiles()
  InstallSamtools()
  InstallRefFlat()
  InstallRefMrnaFile()

  MainSNVClass1(input_file = "lib/data/sample_vcf.txt",
                file_name_in_hla_table = "sample",
                hla_file = "lib/data/sample_hla_table_c1.txt",
                refflat_file  = "lib/refFlat.txt",
                refmrna_file = "lib/refMrna.fa",
                rnaexp_file = "lib/data/sample_rna_exp.txt",
                netMHCpan_dir = "lib/netMHCpan-3.0/netMHCpan",
                nm_id_column = 9
  )

  MainSNVClass2(input_file = "lib/data/sample_vcf.txt",
                file_name_in_hla_table = "sample",
                hla_file = "lib/data/sample_hla_table_c2.txt",
                refflat_file  = "lib/refFlat.txt",
                refmrna_file = "lib/refMrna.fa",
                rnaexp_file = "lib/data/sample_rna_exp.txt",
                netMHCIIpan_dir = "lib/netMHCIIpan-3.1/netMHCIIpan",
                nm_id_column = 9
  )

  MainINDELClass1(input_file = "lib/data/sample_vcf.txt",
                  file_name_in_hla_table = "sample",
                  hla_file = "lib/data/sample_hla_table_c1.txt",
                  refflat_file  = "lib/refFlat.txt",
                  refmrna_file = "lib/refMrna.fa",
                  rnaexp_file = "lib/data/sample_rna_exp.txt",
                  netMHCpan_dir = "lib/netMHCpan-3.0/netMHCpan",
                  nm_id_column = 9
  )

  MainINDELClass2(input_file = "lib/data/sample_vcf.txt",
                  file_name_in_hla_table = "sample",
                  hla_file = "lib/data/sample_hla_table_c2.txt",
                  refflat_file  = "lib/refFlat.txt",
                  refmrna_file = "lib/refMrna.fa",
                  rnaexp_file = "lib/data/sample_rna_exp.txt",
                  netMHCIIpan_dir = "lib/netMHCIIpan-3.1/netMHCIIpan",
                  nm_id_column = 9
  )

  MainMergeSNVClass1(input_dir = "result.sample.NO_job_id_SNV",
                     file_prefix = "NO_job_id_SNV",
                     annotation_file = "lib/data/sample_vcf.txt.NO_job_id_SNV.peptide.txt")

  MainMergeSNVClass2(input_dir = "result.sample.NO_job_id_SNV",
                     file_prefix = "NO_job_id_SNV",
                     annotation_file = "lib/data/sample_vcf.txt.NO_job_id_SNV.peptide.txt")

  MainMergeINDELSVClass1(input_dir = "result.sample.NO_job_id_INDEL",
                         file_prefix = "NO_job_id_INDEL",
                         annotation_file = "lib/data/sample_vcf.txt.NO_JOB_ID_INDEL.peptide.txt")

  MainMergeINDELSVClass2(input_dir = "result.sample.NO_job_id_INDEL",
                         file_prefix = "NO_job_id_INDEL",
                         annotation_file = "lib/data/sample_vcf.txt.NO_JOB_ID_INDEL.peptide.txt")

  MainSVFUSIONClass1(input_file = "lib/data/sample_sv_bnd.tsv",
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

  MainSVFUSIONClass2(input_file = "lib/data/sample_sv_bnd.tsv",
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

  MainMergeINDELSVClass1(input_dir = "result.sample.NO_job_id_SVFusion",
                         file_prefix = "NO_job_id_SVFusion",
                         annotation_file = "lib/data/sample_sv_bnd.tsv.NO_job_id_SVFusion.peptide.txt")

  MainMergeINDELSVClass2(input_dir = "result.sample.NO_job_id_SVFusion",
                         file_prefix = "NO_job_id_SVFusion",
                         annotation_file = "lib/data/sample_sv_bnd.tsv.NO_job_id_SVFusion.peptide.txt")
}
