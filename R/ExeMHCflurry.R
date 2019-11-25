#Execute NetMHCpan
ExemhcflurryClass1<-function(output_peptide_prefix,
                             peptides,
                             hla_types,
                             netMHCpan_dir,
                             peptide_length,
                             export_dir,
                             input_file,
                             job_id){
  print(paste("Executing mhcflurry to", export_dir))
  output_f_header <- ifelse(input_file != "", paste(rev(strsplit(input_file, "/")[[1]])[1], ".", job_id, ".mhcflurry", sep = ""), paste(job_id, ".mhcflurry", sep = ""))
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
      system(paste(MHCflurry,
                   " --mhc-predictor mhcflurry",
                   " --input-fasta-file ", output_f,
                   " --mhc-alleles ", paste("HLA-", gsub("\\*|:","", hla_type), sep = ""),
                   " --mhc-peptide-lengths ", paste(peptide_length, collapse = ","),
                   " --extract-subsequences",
                   " --output-csv ", export_dir, "/", output_f_header, ".HLACLASS1.", COUNT, ".", pep, ".mhcflurry.csv", sep=""))
      COUNT <- COUNT + 1
    }
    if(USETEMP) file.remove(output_f)
  }
}
