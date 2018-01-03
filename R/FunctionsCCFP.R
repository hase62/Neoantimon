CCFP.Calc<-function(cnv_file, output_peptide_txt_file, purity){
  #mutation Rate
  if(ifelse(is.na(cnv_file), FALSE, file.exists(cnv_file))){
    GenerateListForCCFP(output_peptide_txt_file, cnv_file = cnv_file, purity = purity)
    if(file.exists(paste(output_peptide_txt_file, ".cnv.txt", sep=""))){
      calc_ccfp(file_name = paste(output_peptide_txt_file, ".cnv.txt", sep=""))
    }
    GetRatio(output_peptide_txt_file = output_peptide_txt_file,
             output_peptide_txt_cnc_estimate_file = paste(output_peptide_txt_file,".cnv.estimate.txt", sep=""))
  } else {
    GetRatio(output_peptide_txt_file = output_peptide_txt_file,
             output_peptide_txt_cnc_estimate_file = NA)
  }
}

calc_ccfp<-function(file_name, output_ccfp_plot = FALSE){
  prefix <- gsub(".cnv.txt", "", file_name)
  data_ <- scan(file_name, "character", sep="\n")
  if(length(grep("nA|nB|maf|depth|purity", data_)) > 0) data_ <-data_[-1]
  D<-t(sapply(data_, function(x) as.numeric(strsplit(x, "\t")[[1]][-1])))
  id <- sapply(data_, function(x) strsplit(x, "\t")[[1]][1])
  rownames(D)<-id
  cn <- D[,1] + D[,2]
  mac <- D[,3]
  depth <- D[,4]
  purity <- D[,5]
  mac[mac<=1] <- round(mac[mac<=1]*depth[mac<=1])
  #ccf_candidate<-(1:1000)*0.001
  ccf_candidate<-(1:100) * 0.01
  f0 <- purity / (2 * (1 - purity) + cn * purity)
  
  #calculate CCFP
  ccf_posterior<-matrix(nrow = length(f0), ncol = length(ccf_candidate), 0)
  ccf_posterior_cdf<-matrix(nrow = length(f0), ncol = length(ccf_candidate), 0)
  ccf_posterior_cdf5<-rep(0, length(f0)) #lowerBoundmean
  ccf_posterior_cdf50<-rep(0, length(f0)) #mean
  ccf_posterior_cdf95<-rep(0, length(f0)) #upperBound
  colnames(ccf_posterior)<-ccf_candidate
  colnames(ccf_posterior_cdf)<-ccf_candidate
  rownames(ccf_posterior)<-id
  rownames(ccf_posterior_cdf)<-id
  names(ccf_posterior_cdf5)<-id
  names(ccf_posterior_cdf50)<-id
  names(ccf_posterior_cdf95)<-id
  
  for ( i in 1:length(f0) ){
    tmp<-rep(0, length(ccf_candidate))
    f <- f0[i] * ccf_candidate
    f <- ifelse(f > 1, 1, f)
    tmp <- sapply(f, function(x) dbinom(mac[i], depth[i], x))
    tmp2 <- sum(tmp)
    d <- tmp / tmp2
    p <- cumsum(d)
    ccf_posterior[i,] <- d
    ccf_posterior_cdf[i,] <- p
    ccf_posterior_cdf5[i] <- ccf_candidate[min(which(p > 0.05))]
    ccf_posterior_cdf50[i] <- ccf_candidate[min(which(p > 0.5))]
    ccf_posterior_cdf95[i] <- ccf_candidate[min(which(p > 0.95))]
  }
  ccfp_result <- cbind(rownames(D), D[, c(1,2)], ccf_posterior_cdf50, ccf_posterior_cdf5, ccf_posterior_cdf95)
  write.table(x = ccfp_result, file = paste(prefix, ".cnv.estimate.txt", sep=""), 
              row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")
  
  #plot CCFP
  if(output_ccfp_plot ){
    pdf(paste("CCFP", file_name, ".pdf", sep=""))
    plot(ccf_posterior[1,],ylim=c(0,0.3), type="l", ylab="", lwd=0.2)
    for(i in 2:300){
      par(new=T)
      plot(ccf_posterior[i,], ylim=c(0,0.3),type="l", ylab="", lwd=0.2, ann=F)
    }
    
    plot(ccf_posterior_cdf[1,],ylim=c(0,1), type="l", ylab="", lwd=0.2)
    for(i in 2:300){
      par(new=T)
      plot(ccf_posterior_cdf[i,], ylim=c(0,1),type="l", ylab="", lwd=0.2, ann=F)
    }
    
    tmp<-names(sort(ccf_posterior_cdf50))
    plot(ccf_posterior_cdf50[tmp],ylim=c(0,1), pch=20, ylab="",cex = 0.1)
    par(new=T)
    plot(ccf_posterior_cdf5[tmp],ylim=c(0,1),pch=20, ylab="",cex = 0.1, ann=F)
    par(new=T) 
    plot(ccf_posterior_cdf95[tmp],ylim=c(0,1), pch=20, ylab="",cex = 0.1, ann=F)
    dev.off()
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
    #Column + 5
    write.table(cbind(data, matrix(nrow=nrow(data), ncol=5, NA)), 
                output_peptide_txt_file, 
                row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")
  } else if(length(remove)>0){
    #Column + 2
    write.table(cbind(data, cnv[1 + match(data[,1], sapply(cnv[-1,1], function(x) strsplit(x, ";")[[1]][2])),c(2,3)]), 
                  output_peptide_txt_file, 
                  row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")
  } else {
    #Column + 2
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
    ratio_matrix<-matrix(nrow=nrow(data), ncol= 27 - ncol(data), NA)
    write.table(cbind(data, ratio_matrix), output_peptide_txt_file,
                row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")
    return(NULL)
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
    return(NULL)
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
