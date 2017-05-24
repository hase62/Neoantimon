MainINDELClass1<-function(input_file, HLA_file, file_name_in_HLA_table = input_file,
                                       hmdir = getwd(), job_ID = "NO_JOB_ID", RNAseq_file = NA, RNA_bam = NA, 
                                       CNV=NA, Purity = NA,
                                       refDNA = "./../lib_int/GRCh37.fa",
                                       refFlat_file = paste(hmdir,"/../lib_int/refFlat.txt",sep=""), 
                                       refMrna_1 = paste(hmdir,"/../lib_int/refMrna.merge.cut1.fa",sep=""), 
                                       refMrna_3 = paste(hmdir,"/../lib_int/refMrna.merge.cut3.fa",sep=""),
                                       Chr_Column = 1, Mutation_Start_Column = 2, 
                                       Mutation_End_Column = 3, Mutation_Ref_Column = 4, Mutation_Alt_Column = 5, 
                                       NM_ID_Column = 10, Depth_Normal_Column = NA, Depth_Tumor_Column = NA,  
                                       ambiguous_between_exon = 0, ambiguous_codon = 0,
                                       peptide_length = c(8,9,10,11,12,13)){
  
  #Generate FASTA and Mutation Profile
  job_ID = paste(job_ID, "INDEL", sep = "_")
  GenerateIndelSeq(input_file = input_file, hmdir = hmdir, job_ID = job_ID,
                     refFlat_file = refFlat_file, refMrna_1 = refMrna_1, refMrna_3 = refMrna_3,
                     max_peptide_length = max(peptide_length), Chr_Column = Chr_Column, 
                     Mutation_Start_Column = Mutation_Start_Column, 
                     Mutation_End_Column = Mutation_End_Column, 
                     Mutation_Ref_Column = Mutation_Ref_Column, 
                     Mutation_Alt_Column = Mutation_Alt_Column, 
                     NM_ID_Column = NM_ID_Column, 
                     Depth_Normal_Column = Depth_Normal_Column, 
                     Depth_Tumor_Column = Depth_Tumor_Column,
                     ambiguous_between_exon = ambiguous_between_exon, 
                     ambiguous_codon = ambiguous_codon)
  
  output_peptide_txt_file<-paste(input_file,".", job_ID,".peptide.txt",sep="")
  if(!file.exists(output_peptide_txt_file)){
    return(NULL)
  }
    
  #Attach RNAseq Data if Exist, Otherwise set NULL Column
  skip=FALSE
  if(file.exists(RNAseq_file)){
      GenerateListForGetRNASeq(output_peptide_txt_file, width = 120)
      output_file_rna_list<-paste(output_peptide_txt_file, ".list.txt", sep="")
      error<-tryCatch2(system(paste("samtools mpileup -l", output_file_rna_list, "-uf", refDNA, RNA_bam, 
                   ">", paste(output_peptide_txt_file, "list.mp", sep="."))))
      if(error != 0) skip = TRUE
      
      if(!skip) error<-tryCatch2(system(paste("bcftools view -cg", paste(output_peptide_txt_file, "list.mp", sep="."),
                   ">", paste(output_peptide_txt_file, "list.vcf", sep="."))))
      if(error != 0) skip = TRUE
  }
  GetRNAseq_indel(output_peptide_txt_file = output_peptide_txt_file, 
            RNAseq_file = RNAseq_file, 
            output_file_rna_vcf = paste(output_peptide_txt_file, "list.vcf", sep="."))
  
  #Mutation Rate
  if(file.exists(CNV)){
    GenerateListForCCFP(output_peptide_txt_file, CNV = CNV, Purity = Purity)
    system(paste("perl ", gsub("/main","",hmdir), "/lib_MUT/perl/CCFP.pl", " ", 
                 paste(output_peptide_txt_file,".cnv.txt",sep=""), 
           " > ", paste(output_peptide_txt_file,".cnv.estimate.txt",sep=""), sep=""))
    GetRatio(output_peptide_txt_file = output_peptide_txt_file, 
             output_peptide_txt_cnc_estimate_file = paste(output_peptide_txt_file,".cnv.estimate.txt",sep=""))
  }else{
    GetRatio(output_peptide_txt_file = output_peptide_txt_file, 
             output_peptide_txt_cnc_estimate_file = NA)
  }

  #NetMHCpan
  hla<-t(sapply(scan(HLA_file, "character", sep="\n"), function(x) strsplit(x, "\t")[[1]]))
  hla_types<-hla[match(file_name_in_HLA_table, hla[,1]),-1]
  for(pep in c("peptide")){
    COUNT<-1
    for(hla_type in hla_types){
      system(paste("./../lib_int/netMHCpan-3.0/netMHCpan",
                   " -l ", paste(peptide_length, collapse = ","),
                   " -f ", paste(input_file, job_ID, pep,"fasta",sep="."), 
                   " -a HLA-", gsub("\\*","",hla_type), 
                   " > ", input_file, ".", job_ID, ".HLACLASS1.", COUNT, ".", pep, ".txt", 
                   sep=""))
      COUNT <- COUNT + 1
    }
  }
}
                