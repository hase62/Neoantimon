CCFP.Calc<-function(ccfp_dir, cnv_file, output_peptide_txt_file, purity){
  #mutation Rate
  if(ifelse(is.na(ccfp_dir), TRUE, !file.exists(ccfp_dir))) cnv_file<-NA
  if(ifelse(is.na(cnv_file), FALSE, file.exists(cnv_file))){
    GenerateListForCCFP(output_peptide_txt_file, cnv_file = cnv_file, purity = purity)
    print(paste("java -jar ", hmdir, "/", gsub("\\./", "/", ccfp_dir), " ",
                paste(output_peptide_txt_file, ".cnv.txt", sep=""), " > ",
                paste(output_peptide_txt_file, ".cnv.estimate.txt", sep=""), sep=""))
    tryCatch2(system(paste("java -jar ", hmdir, "/", gsub("\\./", "/", ccfp_dir), " ",
                paste(output_peptide_txt_file, ".cnv.txt", sep=""), " > ",
                paste(output_peptide_txt_file, ".cnv.estimate.txt", sep=""), sep="")))
    GetRatio(output_peptide_txt_file = output_peptide_txt_file,
             output_peptide_txt_cnc_estimate_file = paste(output_peptide_txt_file,".cnv.estimate.txt", sep=""))
  } else {
    GetRatio(output_peptide_txt_file = output_peptide_txt_file,
             output_peptide_txt_cnc_estimate_file = NA)
  }
}

GenerateListForCCFP<-function(output_peptide_txt_file, 
                              cnv_file, 
                              purity){
  data<-t(sapply(scan(output_peptide_txt_file, "character", sep="\n"), function(x) strsplit(x, "\t")[[1]]))
  colnames(data)<-NULL
  rownames(data)<-NULL
  if(is.null(nrow(data))){
    data<-rbind(data,data)
  }
  
  data1<-t(sapply(scan(cnv_file, "character", sep="\n"), function(x) strsplit(x, "\t")[[1]]))
  if(length(!is.na(grep("Chr|Pos|R", data1[1,]))) > 0){
    data1<-data1[-1,]
  }
  flag<-FALSE
  if(is.na(purity)){
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
                        data[i,14], data[i,13],
                        purity))
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

GetRatio<-function(output_peptide_txt_file, 
                   output_peptide_txt_cnc_estimate_file){
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
    return()
  }
  
  tmp<-matrix(nrow=nrow(data), ncol=ncol(data_2), NA)
  tmp[match(as.numeric(sapply(data_2[,1], function(x) strsplit(x, ";")[[1]][2])), 1:nrow(data)),]<-data_2
  
  if(nrow(data)==1){
    write.table(t(c(data, tmp[,c(1,2,3,4)])), paste(output_peptide_txt_file,sep=""),
                row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")
  }else{
    write.table(cbind(data, tmp[,c(1,2,3,4)]), paste(output_peptide_txt_file,sep=""),
                row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")
  }
}
