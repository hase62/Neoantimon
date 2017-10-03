GenerateListForCCFP<-function(output_peptide_txt_file, CNV, Purity = NA){
  data<-t(sapply(scan(output_peptide_txt_file, "character", sep="\n"), function(x) strsplit(x, "\t")[[1]]))
  colnames(data)<-NULL
  rownames(data)<-NULL
  if(is.null(nrow(data))){
     data<-rbind(data,data)
  }
  
  data1<-t(sapply(scan(CNV, "character", sep="\n"), function(x) strsplit(x, "\t")[[1]]))
  if(is.na(as.numeric(data1[1,2]))){
   data1<-data1[-1,]
  }
  flag<-FALSE
  if(is.na(Purity)){
   flag<-TRUE
   cnv<-c("","nA","nB","maf","depth","purity")
   remove<-NULL
   cp<-t(sapply(data1[,1], function(x) strsplit(x, ";")[[1]][c(3,4)]))
   for(i in 1:nrow(data)){
    chr<-ifelse(data[i,3]=="X","23", data[i,3])
    chr<-ifelse(data[i,3]=="Y","24",chr)
    dis<-abs(as.numeric(cp[which(chr==cp[,1]), 2]) - as.numeric(data[i,12]))
    j<-1
    data_near<-data1[which(chr==cp[,1])[order(dis)[j]],]
    if(is.na(data_near[6])) next
  
    cnv<-rbind(cnv, c(paste(strsplit(data_near[1],";")[[1]][1], data[i,1], data[i,3], data[i,12],sep=";"), data_near[-1]))
   }
  }
  #read.table(ss, sep=" ")[1,2]
  rownames(data1)<-NULL
  colnames(data1)<-NULL
  
  if(!flag){
   cnv<-c("","nA","nB","maf","depth","purity")
   remove<-NULL
   for(i in 1:nrow(data)){
     chr<-ifelse(data[i,3]=="X","23",data[i,3])
     chr<-ifelse(data[i,3]=="Y","24",chr)
     dis<-abs(as.numeric(data1[which(chr==data1[,2]),3]) - as.numeric(data[i,12]))
     if(length(dis)==0){
        #remove<-c(remove,i)
        next
     }
     j<-1
     while(TRUE){
        j<-j+1
        data_near<-data1[which(chr==data1[,2])[order(dis)[j]],]
        if(data_near[8]!="NA")break
     }
  
     cnv<-rbind(cnv, c(paste(data_near[1], data[i,1], data[i,3], data[i,12],sep=";"),
                       as.numeric(data_near[8]) - as.numeric(data_near[9]),
                       as.numeric(data_near[9]), 
                       data[i,13], data[i,14],
                       Purity))
   }
  }
  #Write List For CCFP
  if(!is.null(nrow(cnv))){
    write.table(cnv, paste(output_peptide_txt_file,".cnv.txt",sep=""), 
       row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")
  }
  #OverWrite peptide.txt
  if(is.null(nrow(cnv))){
    write.table(cbind(data, matrix(nrow=nrow(data), ncol=5, NA)), 
                output_peptide_txt_file, 
                row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")
  }else if(length(remove)>0){
     write.table(cbind(data, cnv[1 + match(data[,1], sapply(cnv[-1,1], function(x) strsplit(x, ";")[[1]][2])),c(2,3)]), 
                 output_peptide_txt_file, 
        row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")
  } else {
     if(nrow(data)==1){
           write.table(t(c(data[1,], cnv[-1,c(2,3)])), output_peptide_txt_file, 
        row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")
     } else {
        write.table(cbind(data, cnv[1 + match(data[,1], sapply(cnv[-1,1], function(x) strsplit(x, ";")[[1]][2])),c(2,3)]), 
                    output_peptide_txt_file, 
           row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")   
     }
  }
}
