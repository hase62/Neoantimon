files<-list.files("NetMHCPan_Results")
files<-files[-grep("_indel", files)]
list.files(paste("NetMHCPan_Results", files[1], sep="/"))

post<-c("MEL","MEL","LUC","LUC","sLUC-SCC","sLUC-AD","BOC","BOC","COC","COC","MAC",
"RCC","RCC","RCC","RCC","HNC","HNC","GAC","GAC","ESC","PRC","PRC","HCC","HCC","HCC",
"HCC","PAC","PAC-END","BLC","THC","BTC","OVC","UTC","CEC","BRC","BRC","BRC",
"sBRC-GBM","sBRC-GBM","sBRC-MED","ML","ML","CLL","AML")

pre<-c("SKCM","MELA","LUSC","LUAD","LUSC","LUAD","BOCA","SARC","COAD","READ","BRCA",
"RECA","KIRC","KIRP","KICH","ORCA","HNSC","STAD","GACA","ESAD","PRAD","EOPC","LINC",
"LIHC","LICA","LIRI","PACA","PAEN","BLCA","THCA","BTCA","OV","UCEC","CESC","GBM",
"LGG","PBCA","GBM","LGG","PBCA","DLBC","MALY","CLLE","LAML")

result<-NULL
co<-0
for(file in files){
   #tumor and HLA
   hit<-grep("HLAtype.txt", list.files(paste("NetMHCPan_Results", file, sep="/")))
   f<-list.files(paste("NetMHCPan_Results", file, sep="/"))[hit]
   hla<-scan(paste("NetMHCPan_Results", file, f, sep="/"), "character", sep="\n")
   #if(length(which(sapply(hla[4:9], nchar)==0))!=0)next
   tumor<-hla[3]
   #count
   hit<-grep("txt.count.txt", list.files(paste("NetMHCPan_Results", file, sep="/")))
   if(length(hit)==0) {
      count<-c(0,0,0,0,0)
   }else{
      f<-list.files(paste("NetMHCPan_Results", file, sep="/"))[hit]
      count<-scan(paste("NetMHCPan_Results", file, f, sep="/"), "character", sep="\n")
   }
   hit<-grep("list.mp", list.files(paste("NetMHCPan_Results", file, sep="/")))
   result<-rbind(result, c(tumor, as.numeric(count[4]), as.numeric(count[5]), length(hit)))
   co<-co+1
}


result<-cbind(result, sapply(result[,1], function(x) strsplit(x, "-")[[1]][1]), 
   post[match(sapply(result[,1], function(x) strsplit(x, "-")[[1]][1]), pre)])
pre_s<-sapply(unique(result[,5]), function(x) sum(as.numeric(result[which(!is.na(match(result[,5], x))),4])))
post_s<-sapply(unique(result[,6]), function(x) sum(as.numeric(result[which(!is.na(match(result[,6], x))),4])))
#write.table(result, "pre_rna.txt", col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\t")
#write.table(result, "post_rna.txt", col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\t")

#d<-scan("result.txt", "character", sep="\n")

