GenerateListForGetRNASeq<-function(output_peptide_txt_file, width = 2){
  data<-t(sapply(scan(output_peptide_txt_file, "character", sep="\n"), 
    function(x) strsplit(x, "\t")[[1]]))
  data<-cbind(data[,3],as.numeric(data[,12]) - width, as.numeric(data[,12]) + width)
  write.table(data, paste(output_peptide_txt_file, ".list.txt", sep=""), row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")
}
