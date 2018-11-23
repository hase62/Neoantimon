#Execute NetMHCpan
ExeNetMHCpanClass1<-function(output_peptide_prefix,
                             peptides,
                             hla_types,
                             netMHCpan_dir,
                             peptide_length,
                             export_dir,
                             input_file,
                             job_id){
  print(paste("Executing netMHCpan to", export_dir))
  for(pep in peptides){
    COUNT<-1
    output_f <- paste(output_peptide_prefix, pep, "fasta",sep=".")
    USETEMP <- FALSE
    if(nchar(output_f) > 230) {
      output_f_new <- paste("temp.Neoantimon.", runif(1) * 1000000, "txt", sep = "")
      file.copy(from = output_f, to = output_f_new)
      output_f <- output_f_new
      USETEMP <- TRUE
    }
    for(hla_type in hla_types){
      paste("Calculating", pep, hla_type)
      system(paste(netMHCpan_dir,
                   " -BA ",
                   " -l ", paste(peptide_length, collapse = ","),
                   " -f ", output_f,
                   " -a HLA-", gsub("\\*","",hla_type),
                   " > ", export_dir, "/", rev(strsplit(input_file, "/")[[1]])[1], ".", job_id, ".HLACLASS1.", COUNT, ".", pep, ".txt", sep=""))
      COUNT <- COUNT + 1
    }
    if(USETEMP) file.remove(output_f)
  }
}

ExeNetMHCpanClass2<-function(output_peptide_prefix,
                               peptides,
                               hla_types,
                               netMHCIIpan_dir,
                               peptide_length,
                               export_dir,
                               input_file,
                               job_id){
  print(paste("Executing netMHCpan to", export_dir))
  for(pep in peptides){
    COUNT<-1
    output_f <- paste(output_peptide_prefix, pep, "fasta",sep=".")
    USETEMP <- FALSE
    if(nchar(output_f) > 230) {
      output_f_new <- paste("temp.Neoantimon.", runif(1) * 1000000, "txt", sep = "")
      file.copy(from = output_f, to = output_f_new)
      output_f <- output_f_new
      USETEMP <- TRUE
    }
    for(hla_type in hla_types){
      if(length(grep("DRB1", hla_type))==1) {
        system(paste(netMHCIIpan_dir,
                     " -length ", paste(peptide_length, collapse = ","),
                     " -f ", output_f,
                     " -a ", gsub("\\*","_", gsub("\\:","",hla_type)),
                     " > ", export_dir, "/", rev(strsplit(input_file, "/")[[1]])[1], ".", job_id, ".HLACLASS2.", COUNT, ".", pep, ".txt", sep=""))
        COUNT <- COUNT + 1
      }

      if(length(grep("DPA1", hla_type))==1) {
        for(hla2 in hla_types[grep("DPB1", hla_types)]){
          system(paste(netMHCIIpan_dir,
                       " -length ", paste(peptide_length, collapse = ","),
                       " -f ", output_f,
                       " -choose -cha ", gsub("\\*|\\:","", hla_type),
                       " -choose -chb ", gsub("\\*|\\:","", hla2),
                       " > ", export_dir, "/", rev(strsplit(input_file, "/")[[1]])[1], ".", job_id, ".HLACLASS2.", COUNT, ".", pep, ".txt", sep=""))
          COUNT <- COUNT + 1
        }
      }

      if(length(grep("DQA1", hla_type))==1) {
        for(hla2 in hla_types[grep("DQB1", hla_types)]){
          system(paste(netMHCIIpan_dir,
                       " -length ", paste(peptide_length, collapse = ","),
                       " -f ", output_f,
                       " -choose -cha ", gsub("\\*|\\:","", hla_type),
                       " -choose -chb ", gsub("\\*|\\:","", hla2),
                       " > ", export_dir, "/", rev(strsplit(input_file, "/")[[1]])[1], ".", job_id, ".HLACLASS2.", COUNT, ".", pep, ".txt", sep=""))
          COUNT <- COUNT + 1
        }
      }
    }
    if(USETEMP) file.remove(output_f)
  }
}
