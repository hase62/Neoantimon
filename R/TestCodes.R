#Initial Setting
file_name_in_hla_table = input_file
hmdir = getwd() 
refflat_file = paste(hmdir, "/lib/refFlat.txt", sep="")
refmrna_file = paste(hmdir, "/lib/refMrna.merge.fa", sep="")
job_id = "NO_job_id" 
rnaexp_file = NA 
rnabam_file = NA
cnv_file=NA 
ccfp_dir = paste(hmdir, "lib/ccfp.jar", sep="") 
purity = 1
netMHCpan_dir = paste(hmdir, "/lib/netMHCpan-3.0/netMHCIIpan", sep="")
netMHCIIpan_dir = paste(hmdir, "/lib/netMHCIIpan-3.1/netMHCIIpan", sep="")
refdna_file = paste(hmdir, "/lib/GRCh37.fa", sep="")
samtools_dir = NA
bcftools_dir = NA
chr_column = NA 
mutation_start_column = NA
mutation_end_column = NA 
mutation_ref_column = NA 
mutation_alt_column = NA
nm_id_column = NA 
depth_normal_column = NA 
depth_tumor_column = NA
ambiguous_between_exon = 0 
ambiguous_codon = 0
peptide_length = c(8, 9, 10, 11, 12, 13)

hmdir = getwd()
input_file = "lib/data/sample_vcf.txt"
job_id  = "NO_JOB_ID"
file_name_in_hla_table = "sample"
hla_file = "lib/data/sample_hla_table_c1.txt"
refflat_file  = "lib/refFlat.txt"
refmrna_file = "lib/refMrna.merge.fa"
rnaexp_file = "lib/data/sample_rna_exp.txt"
cnv_file = "lib/data/sample_copynum.txt"
purity = 0.8
samtools_dir = "lib/samtools-0.1.19/samtools"
bcftools_dir = "lib/samtools-0.1.19/bcftools/bcftools"
ccfp_dir = "lib/ccfp.jar"
netMHCpan_dir = "lib/netMHCpan-3.0/netMHCpan"

files<-list.files()
for(file in files){
  source(file)
}

MainSNVClass1(input_file = "lib/data/sample_vcf.txt",
              file_name_in_hla_table = "sample",
              hla_file = "lib/data/sample_hla_table_c1.txt",
              refflat_file  = "lib/refFlat.txt",
              refmrna_file = "lib/refMrna.merge.fa",
              rnaexp_file = "lib/data/sample_rna_exp.txt",
              cnv_file = "lib/data/sample_copynum.txt",
              ccfp_dir = "lib/ccfp.jar",
              netMHCpan_dir = "lib/netMHCpan-3.0/netMHCpan",
              refdna_file = "lib/GRCh37.fa"
              )

MainINDELClass1(input_file = "lib/data/sample_vcf.txt",
              file_name_in_hla_table = "sample",
              hla_file = "lib/data/sample_hla_table_c1.txt",
              refflat_file  = "lib/refFlat.txt",
              refmrna_file = "lib/refMrna.merge.fa",
              rnaexp_file = "lib/data/sample_rna_exp.txt",
              cnv_file = "lib/data/sample_copynum.txt",
              ccfp_dir = "lib/ccfp.jar",
              netMHCpan_dir = "lib/netMHCpan-3.0/netMHCpan",
              refdna_file = "lib/GRCh37.fa"
              )

hla_file = "lib/data/sample_hla_table_c2.txt"
netMHCIIpan_dir = "lib/netMHCIIpan-3.1/netMHCIIpan"
peptide_length = c(15)

MainSNVClass2(input_file = "lib/data/sample_vcf.txt",
              file_name_in_hla_table = "sample",
              hla_file = "lib/data/sample_hla_table_c2.txt",
              refflat_file  = "lib/refFlat.txt",
              refmrna_file = "lib/refMrna.merge.fa",
              rnaexp_file = "lib/data/sample_rna_exp.txt",
              cnv_file = "lib/data/sample_copynum.txt",
              ccfp_dir = "lib/ccfp.jar",
              netMHCIIpan_dir = "lib/netMHCIIpan-3.1/netMHCIIpan"
)

MainINDELClass2(input_file = "lib/data/sample_vcf.txt",
              file_name_in_hla_table = "sample",
              hla_file = "lib/data/sample_hla_table_c2.txt",
              refflat_file  = "lib/refFlat.txt",
              refmrna_file = "lib/refMrna.merge.fa",
              rnaexp_file = "lib/data/sample_rna_exp.txt",
              cnv_file = "lib/data/sample_copynum.txt",
              ccfp_dir = "lib/ccfp.jar",
              netMHCIIpan_dir = "lib/netMHCIIpan-3.1/netMHCIIpan"
)


input_dir <- "result.sample.NO_job_id_SNV"
file_prefix <- "NO_job_id_SNV"
annotation_file<-"lib/data/sample_vcf.txt.NO_job_id_SNV.peptide.txt"

MainMergeSNVClass1(input_dir = input_dir,
                   file_prefix = file_prefix,
                   annotation_file = annotation_file)

MainMergeSNVClass2(input_dir = input_dir,
                   file_prefix = file_prefix,
                   annotation_file = annotation_file)

input_dir <- "result.sample.NO_job_id_INDEL/"
file_prefix <- "NO_job_id_INDEL"
annotation_file<-"lib/data/sample_vcf.txt.NO_JOB_ID_INDEL.peptide.txt"

MainMergeINDELClass1(input_dir = input_dir,
                   file_prefix = file_prefix,
                   annotation_file = annotation_file)

MainMergeINDELClass2(input_dir = input_dir,
                   file_prefix = file_prefix,
                   annotation_file = annotation_file)




