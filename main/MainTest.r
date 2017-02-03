#Class1
source("MainClass1.r")
MainNeoantigenidentification(hmdir = getwd(), 
                             input_file="./../lib_sample/sample_annovar.txt", 
                             job_ID = "NO_JOB_ID", 
                             file_name_in_HLA_table = "096b4f32-10c1-4737-a0dd-cae04c54ee33", 
                             HLA_file = "./../lib_sample/hla_table.txt", 
                             RNAseq_file = "./../lib_sample/RNAseq.txt", 
                             RNA_bam="./../lib_sample/RNAbam.bam", 
                             Genomon = FALSE, 
                             CNV="./../lib_sample/Copy.txt", 
                             Purity = 0.719526227140365,
                             refDNA = "./../lib_int/GRCh37.fa")

