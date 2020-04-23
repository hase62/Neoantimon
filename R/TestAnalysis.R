#'Execute Sample Analysis
#'
#'@return void
#'
#'@export
TestAnalysis<-function(){
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

}
