files<-list.files(commandArgs(TRUE)[1])

post<-c("MEL","MEL","LUC","LUC","sLUC-SCC","sLUC-AD","BOC","BOC","COC","COC","MAC",
"RCC","RCC","RCC","RCC","HNC","HNC","GAC","GAC","ESC","PRC","PRC","HCC","HCC","HCC",
"HCC","PAC","PAC-END","BLC","THC","BTC","OVC","UTC","CEC","BRC","BRC","BRC",
"sBRC-GBM","sBRC-GBM","sBRC-MED","ML","ML","CLL","AML")

pre<-c("SKCM","MELA","LUSC","LUAD","LUSC","LUAD","BOCA","SARC","COAD","READ","BRCA",
"RECA","KIRC","KIRP","KICH","ORCA","HNSC","STAD","GACA","ESAD","PRAD","EOPC","LINC",
"LIHC","LICA","LIRI","PACA","PAEN","BLCA","THCA","BTCA","OV","UCEC","CESC","GBM",
"LGG","PBCA","GBM","LGG","PBCA","DLBC","MALY","CLLE","LAML")

up<-unique(pre)
cinfo<-c("", "Gene ID", "Chr", "NM_ID", "Change", "ref", "alt", "Prob", "Mutation Prob.",
                "Exon Start", "Exon End", "Mutation Position", "Depth", "TumorDepth",
                "Peptide Normal", "Peptide Mutation", "DNA_Normal", "DNA_Mut",
                "TotalRNA", "TumorRNARatio", "TumorRNA", "nA", "nB", "Checker",
                "MutRatio", "MutRatio Min", "MutRatio Max")

result<-list(NULL)
for(i in 1:length(up)){
   result<-c(result, list(NULL))
   
}
data2<-NULL
miss<-NULL
for(file in files){
   #tumor and HLA
   fs<-list.files(paste(commandArgs(TRUE)[1], file, sep="/"))
   hit<-grep("HLAtype.txt", fs)
   f<-fs[hit]
   if(length(f)==0) next
   hla<-scan(paste(commandArgs(TRUE)[1], file, f, sep="/"), "character", sep="\n")[3]
   #if(is.na(match(strsplit(hla,"-")[[1]][1],lung)))next 
   if(is.na(fs[grep("output.tsv.peptide.indel.txt",fs)][1])) {
      next
   }
   print(hla)
   data<-t(sapply(scan(paste(commandArgs(TRUE)[1], file, sep="/", fs[grep("output.tsv.peptide.indel.txt",fs)][1]),"character",sep="\n"), 
      function(x) strsplit(x,"\t")[[1]]))
   if(ncol(data)!=27) {
      miss<-c(miss,file)
      next
   }
   colnames(data)<-NULL
   rownames(data)<-NULL
   data<-data[nchar(data[,15])!=1,]
   if(is.null(dim(data))){data<-t(as.matrix(data))}
   data<-cbind(f,strsplit(hla,"-")[[1]][1], data,c(nchar(data[,16]) - nchar(data[,15])))
   data2<-rbind(data2, data)
}
print(miss)

###Remove Tumor
unq<-sort(as.numeric(unique(data2[,30])))
depth<-sapply(data2[,22], function(x) strsplit(x, "/")[[1]][2])
len_m<-sapply(data2[,19], nchar)
len_n<-sapply(data2[,20], nchar)
mRNAl<-sapply(1:nrow(data2), function(x) strsplit(data2[x,19],"")[[1]] == strsplit(data2[x,20],"")[[1]])
len<-unlist(lapply(1:length(mRNAl), function(x) nchar(data2[x,20]) - which(!mRNAl[[x]])[1]))
data2<-cbind(data2,len_m,len)


###Plot
for(s in c("alpha", "beta")){
 pdf(paste("FullPlot_", s, "_tumor.pdf", sep=""))
 for(n in c(-1, 0, 1, 5, 10, 30)){
  for(i in 1:length(unique(post))){
   for(cut in c(50, 100, 150, 200)){
    for(j in c(23)){
     if(j==21) {tag<-"Total"
     } else {tag<-"Tumor"}
     hit<-which(!is.na(match(data2[,2], pre[which(!is.na(match(post,unique(post)[i])))])))
     ##FPKM
     hit2<-hit[data2[hit, j]!="NA"]
     if(length(hit2)==0) next
     if(s=="beta") hit2<-hit2[as.numeric(depth[hit2])>5]
     if(length(hit2)==0) next
     hit2<-hit2[which(as.numeric(data2[hit2,j])>n)]
    
     result<-cbind(log(as.numeric(data2[hit2,30])+1), log(as.numeric(data2[hit2,j])+1))
     result2<-cbind(log(as.numeric(data2[hit2,31])+1), log(as.numeric(data2[hit2,j])+1))
     result3<-cbind(log(as.numeric(data2[hit2,32])+1), log(as.numeric(data2[hit2,j])+1))
     if(length(result[exp(result[,1])<=cut, 1])==0) next
     if(length(result[exp(result[,1])>cut, 1])==0) next
     test<-wilcox.test(result[exp(result[,1])<=cut, 2], result[exp(result[,1])>cut, 2])
     if(cut==50) plot(result2, xlab="mRNA Full Length", ylab="Exp", main=paste(unique(post)[i], tag, n, "mRNA Full Length"))
     if(cut==50) plot(result3, xlab="mRNA Red. Length", ylab="Exp", main=paste(unique(post)[i], tag, n, "mRNA Red Length"))
     plot(result, xlab="length", ylab="Exp", main=paste(unique(post)[i], tag, n, "p=", test$p.value))
     abline(v = log(cut+1))
     if(n>-1) abline(h = log(n+1))
     if(n > -1){
      boxplot(result[exp(result[,1])<=cut, 2], result[exp(result[,1])>cut, 2], 
	   main=paste(unique(post)[i], tag, n, "p=", test$p.value))
     }
    }
   }
  }
 }
 dev.off()
}

pdf("IndelLength_post_.pdf")
th<-200
depth_th<-5
for(i in 1:length(unique(post))){
for(exp in c(0,1,2,3,4,5)){
  for(j in c(23)){
   if(j==21) {tag<-"Total"
   } else {tag<-"Tumor"}
   hit<-which(!is.na(match(data2[,2], pre[which(!is.na(match(post,unique(post)[i])))])))
   hit<-hit[as.numeric(data2[hit, 30])<th]
   ##All
   tmp<-as.numeric(data2[hit, 30])
   s1<-length(tmp)
   count<-NULL
   for(u in unique(tmp)){count<-rbind(count, c(u, length(which(u==tmp))))}
   count<-count[order(count[,1]),]
   mat<-matrix(nrow=th,ncol=2,0)
   mat[,1]<-seq(1:th)-1
   mat[match(count[,1], mat[,1]), 2]<-count[,2]
   mat<-cbind(mat, mat[,2]/sum(mat[,2]),(61/64)^(mat[,1]) * 3/64)
   mat<-cbind(mat, log(mat[,3]))
   ##FPKM
   hit2<-hit[data2[hit, j]!="NA"]
   tmp<-as.numeric(data2[hit2[as.numeric(data2[hit2,j])>exp & as.numeric(depth[hit2]) > depth_th], 30])
   s2<-length(tmp)
   if(length(tmp)>1){
      count<-NULL
      for(u in unique(tmp)){count<-rbind(count, c(u, length(which(u==tmp))))}
      count<-count[order(count[,1]),]
      mat2<-matrix(nrow=th,ncol=2,0)
      mat2[,1]<-seq(1:th)-1
      mat2[match(count[,1], mat2[,1]), 2]<-count[,2]
      mat2<-cbind(mat2, mat2[,2]/sum(mat2[,2]),(61/64)^(mat2[,1]) * 3/64)
      mat2<-cbind(mat2, log(mat2[,3]))
   } else {mat2<-matrix(nrow=th,ncol=5,0)}
   ##FPKM
   tmp<-as.numeric(data2[hit2[as.numeric(data2[hit2,j])<=exp & as.numeric(data2[hit2,j])>=0 & as.numeric(depth[hit2]) > depth_th], 30])
   s3<-length(tmp)
   if(length(tmp)>1){
      count<-NULL
      for(u in unique(tmp)){count<-rbind(count, c(u, length(which(u==tmp))))}
      if(nrow(count)==1) {count<-t(count[order(count[,1]),])
      } else {count<-count[order(count[,1]),]}
      mat3<-matrix(nrow=th,ncol=2,0)
      mat3[,1]<-seq(1:th)-1
      mat3[match(count[,1], mat3[,1]), 2]<-count[,2]
      mat3<-cbind(mat3, mat3[,2]/sum(mat3[,2]),(61/64)^(mat3[,1]) * 3/64)
      mat3<-cbind(mat3, log(mat3[,3]))
   } else {mat3<-matrix(nrow=th,ncol=5,0)}
   matplot(cbind(mat[,c(4,3)],mat2[,3],mat3[,3]), type="l",
      main = paste(unique(post)[i], tag, ": th=", exp), xlim=c(0,th),lwd=1,lty=c(1,1,1,1),ylab="Rate")
   legend("topright",legend=c("Theory", paste("All:n=", s1), paste("x>", exp, ":n=", s2), paste("0<x<=", exp, ":n=", s3)), lwd=1, lty=1, col=1:4)
   matplot(cbind(log(mat[,4]), mat[,5], mat2[,5], mat3[,5]), type="l", main = paste(unique(post)[i], tag, ": th=", exp),
      xlim=c(0,th),lwd=1,lty=c(1,1,1,1),ylab="Rate")
   legend("topright",legend=c("Theory", paste("All:n=", s1), paste("x>", exp, ":n=", s2), paste("0<x<=", exp, ":n=", s3)), lwd=1, lty=1, col=1:4)
  }
}
}
dev.off()

