source("Functions/MainNeoantimon.r")

##
#lib_sample includes all sample files.
##
MainSNVClass1(hmdir = getwd(),
              input_file = "./../lib_sample/sample_annovar.txt", 
              job_ID = "NO_JOB_ID", 
              file_name_in_HLA_table = "sample", 
              HLA_file = "./../lib_sample/hla_table.txt", 
              RNAseq_file = "./../lib_sample/RNAseq.txt", 
              RNA_bam="./../lib_sample/RNAbam.bam", 
              CNV="./../lib_sample/Copy.txt", 
              Purity = 0.8,
              refDNA = "./../lib_int/GRCh38.fa")

MainMergeClass1(hmdir = getwd(), 
                input_dir = "./../lib_sample", 
                input_file_prefix = "sample_annovar", 
                Tumor_RNA_BASED_ON_DNA = TRUE)

MainSNVClass1(hmdir = getwd(),
              input_file = "./../lib_sample/sample_genomon.txt", 
              job_ID = "NO_JOB_ID", 
              file_name_in_HLA_table = "sample",
              HLA_file = "./../lib_sample/hla_table.txt", 
              RNAseq_file = "./../lib_sample/RNAseq.txt", 
              RNA_bam="./../lib_sample/RNAbam.bam", 
              CNV="./../lib_sample/Copy.txt", 
              Purity = 0.8,
              refDNA = "./../lib_int/GRCh38.fa",
              Depth_Normal_Column = 53, Depth_Tumor_Column = 51)

MainMergeClass1(hmdir = getwd(), 
                input_dir = "./../lib_sample", 
                input_file_prefix = "sample_genomon", 
                Tumor_RNA_BASED_ON_DNA = TRUE)

MainSNVClass2(hmdir = getwd(), 
              input_file="./../lib_sample/sample_annovar.txt", 
              job_ID = "NO_JOB_ID", 
              file_name_in_HLA_table = "sample", 
              HLA_file = "./../lib_sample/hla_table2.txt", 
              RNAseq_file = "./../lib_sample/RNAseq.txt", 
              RNA_bam="./../lib_sample/RNAbam.bam", 
              CNV="./../lib_sample/Copy.txt", 
              Purity = 0.8,
              refDNA = "./../lib_int/GRCh38.fa")
MainMergeClass2(hmdir = getwd(), 
                input_dir = "./../lib_sample", 
                input_file_prefix = "sample_annovar", 
                Tumor_RNA_BASED_ON_DNA = TRUE)

MainSNVClass2(hmdir = getwd(), 
              input_file="./../lib_sample/sample_genomon.txt", 
              job_ID = "NO_JOB_ID", 
              file_name_in_HLA_table = "sample",
              HLA_file = "./../lib_sample/hla_table2.txt", 
              RNAseq_file = "./../lib_sample/RNAseq.txt", 
              RNA_bam="./../lib_sample/RNAbam.bam", 
              CNV="./../lib_sample/Copy.txt", 
              Purity = 0.8,
              refDNA = "./../lib_int/GRCh38.fa",
              Depth_Normal_Column = 53, Depth_Tumor_Column = 51)

MainMergeClass2(hmdir = getwd(), 
                input_dir = "./../lib_sample", 
                input_file_prefix = "sample_genomon", 
                Tumor_RNA_BASED_ON_DNA = TRUE)

MainINDELClass1(hmdir = getwd(),
              input_file = "./../lib_sample/sample_genomon.txt", 
              job_ID = "NO_JOB_ID", 
              file_name_in_HLA_table = "sample", 
              HLA_file = "./../lib_sample/hla_table.txt", 
              RNAseq_file = "./../lib_sample/RNAseq.txt", 
              RNA_bam="./../lib_sample/RNAbam.bam", 
              CNV="./../lib_sample/Copy.txt", 
              Purity = 0.8,
              refDNA = "./../lib_int/GRCh37.fa",
              Depth_Normal_Column = 53, Depth_Tumor_Column = 51)

MainMergeClass1(hmdir = getwd(), 
                input_dir = "./../lib_sample", 
                input_file_prefix = "sample_genomon", 
                Tumor_RNA_BASED_ON_DNA = TRUE,
                INDEL = TRUE)

MainINDELClass2(hmdir = getwd(),
                input_file = "./../lib_sample/sample_genomon.txt", 
                job_ID = "NO_JOB_ID", 
                file_name_in_HLA_table = "sample",
                HLA_file = "./../lib_sample/hla_table2.txt", 
                RNAseq_file = "./../lib_sample/RNAseq.txt", 
                RNA_bam="./../lib_sample/RNAbam.bam", 
                CNV="./../lib_sample/Copy.txt", 
                Purity = 0.8,
                refDNA = "./../lib_int/GRCh37.fa",
                Depth_Normal_Column = 53, Depth_Tumor_Column = 51)

MainMergeClass2(hmdir = getwd(), 
                input_dir = "./../lib_sample", 
                input_file_prefix = "sample_genomon", 
                Tumor_RNA_BASED_ON_DNA = TRUE, 
                INDEL = TRUE)

