GetRatio<-function(output_peptide_txt_file, output_peptide_txt_cnc_estimate_file = NA){
  data<-t(sapply(scan(output_peptide_txt_file, "character", sep="\n"), function(x) strsplit(x, "\t")[[1]]))
  colnames(data)<-NULL
  rownames(data)<-NULL
  if(ifelse(is.na(output_peptide_txt_cnc_estimate_file), TRUE, !file.exists(output_peptide_txt_cnc_estimate_file))){
     ratio_matrix<-matrix(nrow=nrow(data), ncol=6, NA)
     write.table(cbind(data, ratio_matrix), output_peptide_txt_file,
        row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")
     return()
  }

  #name, mean, min, max
  data_2<-t(sapply(scan(output_peptide_txt_cnc_estimate_file, "character", sep="\n", skip=1),
                   function(x) strsplit(x, "\t")[[1]]))
  colnames(data_2)<-NULL
  rownames(data_2)<-NULL
  if(length(data_2)==0){
     ratio_matrix<-matrix(nrow=nrow(data), ncol= 27 - ncol(data), NA)
     write.table(cbind(data, ratio_matrix), output_peptide_txt_file,
        row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")
     q("no")
  }

  tmp<-matrix(nrow=nrow(data), ncol=ncol(data_2), -1)
  tmp[match(as.numeric(sapply(data_2[,1], function(x) strsplit(x, ";")[[1]][2])), 1:nrow(data)),]<-data_2

  if(nrow(data)==1){
     write.table(t(c(data, tmp[,c(1,2,3,4)])), paste(output_peptide_txt_file,sep=""),
        row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")
  }else{
     write.table(cbind(data, tmp[,c(1,2,3,4)]), paste(output_peptide_txt_file,sep=""),
        row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")
  }
}
