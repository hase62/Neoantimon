GenerateListForGetRNASeq<-function(output_peptide_txt_file){
  data<-t(sapply(scan(output_peptide_txt_file, "character", sep="\n"), 
  		function(x) strsplit(x, "\t")[[1]]))
  data<-cbind(data[,3],as.numeric(data[,12]) - 2,as.numeric(data[,12]) + 2)
  write.table(data, paste(output_peptide_txt_file, ".list.txt", sep=""), row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")
}