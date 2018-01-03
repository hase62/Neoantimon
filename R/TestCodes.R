#Initial Setting
samtools_dir = NA
bcftools_dir = NA
chr_column = NA 
mutation_start_column = NA
mutation_end_column = NA 
mutation_ref_column = NA 
mutation_alt_column = NA
nm_id_column = NA 
mutation_alt_bnd_column = 18
gene_symbol_column = 7
mate_id_column<-16
depth_normal_column = NA 
depth_tumor_column = NA 
ambiguous_between_exon = 0 
ambiguous_codon = 0 
peptide_length = c(8, 9, 10, 11, 12, 13) 

hmdir = getwd()
input_file = "lib/data/sample_vcf.txt"
file_name_in_hla_table = "sample"
job_id  = "NO_JOB_ID"
file_name_in_hla_table = "sample"
hla_file = "lib/data/sample_hla_table_c1.txt"
refflat_file  = "lib/refFlat.txt"
refmrna_file = "lib/refMrna.fa"
refdna_file = paste(hmdir, "/lib/GRCh37.fa", sep="")

rnabam_file = NA
purity = 1

rnaexp_file = NA 
rnaexp_file = "lib/data/sample_rna_exp2.txt"
cnv_file=NA 
cnv_file = "lib/data/sample_copynum2.txt"
purity = 0.8
samtools_dir = "lib/samtools-0.1.19/samtools"
bcftools_dir = "lib/samtools-0.1.19/bcftools/bcftools"
netMHCpan_dir = paste(hmdir, "/lib/netMHCpan-3.0/netMHCpan", sep="")
netMHCIIpan_dir = paste(hmdir, "/lib/netMHCIIpan-3.1/netMHCIIpan", sep="")

files<-list.files()
files<-files[grep(".R", files)]
files<-files[grep("TestCodes", files, invert = TRUE)]
for(file in files){
  source(file)
}

MainSNVClass1(input_file = "lib/data/sample_vcf.txt",
              file_name_in_hla_table = "sample",
              hla_file = "lib/data/sample_hla_table_c1.txt",
              refflat_file  = "lib/refFlat.txt",
              refmrna_file = "lib/refMrna.fa",
              rnaexp_file = "lib/data/sample_rna_exp2.txt",
              cnv_file = "lib/data/sample_copynum2.txt",
              netMHCpan_dir = "lib/netMHCpan-3.0/netMHCpan",
              refdna_file = "lib/GRCh37.fa"
              )

MainSNVClass1(input_file = "lib/data/sample_vcf2.txt",
              file_name_in_hla_table = "sample2",
              hla_file = "lib/data/sample_hla_table_c1.txt",
              refflat_file  = "lib/refFlat.txt",
              refmrna_file = "lib/refMrna.fa",
              rnaexp_file = "lib/data/sample_rna_exp2.txt",
              cnv_file = "lib/data/sample_copynum2.txt",
              netMHCpan_dir = "lib/netMHCpan-3.0/netMHCpan",
              refdna_file = "lib/GRCh37.fa",
              rnabam_file = "lib/RNAbam.bam",
              nm_id_column = 10,
              depth_normal_column = 11,
              depth_tumor_column = 12
)

MainINDELClass1(input_file = "lib/data/sample_vcf.txt",
              file_name_in_hla_table = "sample",
              hla_file = "lib/data/sample_hla_table_c1.txt",
              refflat_file  = "lib/refFlat.txt",
              refmrna_file = "lib/refMrna.fa",
              rnaexp_file = "lib/data/sample_rna_exp2.txt",
              cnv_file = "lib/data/sample_copynum2.txt",
              netMHCpan_dir = "lib/netMHCpan-3.0/netMHCpan",
              refdna_file = "lib/GRCh37.fa"
              )

MainSVFUSIONClass1(input_file = "lib/data/sample_sv_bnd.tsv",
              file_name_in_hla_table = "sample",
              hla_file = "lib/data/sample_hla_table_c1.txt",
              refflat_file  = "lib/refFlat.txt",
              refmrna_file = "lib/refMrna.fa",
              rnaexp_file = "lib/data/sample_rna_exp2.txt",
              cnv_file = "lib/data/sample_copynum2.txt",
              netMHCpan_dir = "lib/netMHCpan-3.0/netMHCpan",
              refdna_file = "lib/GRCh37.fa",
              mutation_alt_bnd_column = 18,
              gene_symbol_column = 7,
              mate_id_column = 16
)

hla_file = "lib/data/sample_hla_table_c2.txt"
netMHCIIpan_dir = "lib/netMHCIIpan-3.1/netMHCIIpan"
peptide_length = c(15)

MainSNVClass2(input_file = "lib/data/sample_vcf.txt",
              file_name_in_hla_table = "sample",
              hla_file = "lib/data/sample_hla_table_c2.txt",
              refflat_file  = "lib/refFlat.txt",
              refmrna_file = "lib/refMrna.fa",
              rnaexp_file = "lib/data/sample_rna_exp2.txt",
              cnv_file = "lib/data/sample_copynum2.txt",
              netMHCIIpan_dir = "lib/netMHCIIpan-3.1/netMHCIIpan"
)

MainINDELClass2(input_file = "lib/data/sample_vcf.txt",
              file_name_in_hla_table = "sample",
              hla_file = "lib/data/sample_hla_table_c2.txt",
              refflat_file  = "lib/refFlat.txt",
              refmrna_file = "lib/refMrna.fa",
              rnaexp_file = "lib/data/sample_rna_exp2.txt",
              cnv_file = "lib/data/sample_copynum2.txt",
              netMHCIIpan_dir = "lib/netMHCIIpan-3.1/netMHCIIpan"
)

MainSVFUSIONClass2(input_file = "lib/data/sample_sv_bnd.tsv",
                   file_name_in_hla_table = "sample",
                   hla_file = "lib/data/sample_hla_table_c2.txt",
                   refflat_file  = "lib/refFlat.txt",
                   refmrna_file = "lib/refMrna.fa",
                   rnaexp_file = "lib/data/sample_rna_exp2.txt",
                   cnv_file = "lib/data/sample_copynum2.txt",
                   netMHCIIpan_dir = "lib/netMHCIIpan-3.1/netMHCIIpan",
                   refdna_file = "lib/GRCh37.fa",
                   mutation_alt_bnd_column = 18,
                   gene_symbol_column = 7,
                   mate_id_column = 16
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

MainMergeINDELSVClass1(input_dir = input_dir,
                   file_prefix = file_prefix,
                   annotation_file = annotation_file)

MainMergeINDELSVClass2(input_dir = input_dir,
                   file_prefix = file_prefix,
                   annotation_file = annotation_file)

input_dir <- "result.sample.NO_job_id_SVFusion"
file_prefix <- "SVFusion"
annotation_file<-"lib/data/sample_sv_bnd.tsv.NO_job_id_SVFusion.peptide.txt"

MainMergeINDELSVClass1(input_dir = input_dir,
                     file_prefix = file_prefix,
                     annotation_file = annotation_file)

MainMergeINDELSVClass2(input_dir = input_dir,
                     file_prefix = file_prefix,
                     annotation_file = annotation_file)



