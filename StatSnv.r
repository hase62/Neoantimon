root<-"Output.Snv.1"
files<-list.files(root)

post<-c("MEL","MEL","LUC","LUC","sLUC-SCC","sLUC-AD","BOC","BOC","COC","COC","MAC",
"RCC","RCC","RCC","RCC","HNC","HNC","GAC","GAC","ESC","PRC","PRC","HCC","HCC","HCC",
"HCC","PAC","PAC-END","BLC","THC","BTC","OVC","UTC","CEC","BRC","BRC","BRC",
"sBRC-GBM","sBRC-GBM","sBRC-MED","ML","ML","CLL","AML")

pre<-c("SKCM","MELA","LUSC","LUAD","LUSC","LUAD","BOCA","SARC","COAD","READ","BRCA",
"RECA","KIRC","KIRP","KICH","ORCA","HNSC","STAD","GACA","ESAD","PRAD","EOPC","LINC",
"LIHC","LICA","LIRI","PACA","PAEN","BLCA","THCA","BTCA","OV","UCEC","CESC","GBM",
"LGG","PBCA","GBM","LGG","PBCA","DLBC","MALY","CLLE","LAML")

cinfo<-c("", "Gene ID", "Chr", "NM_ID", "Change", "ref", "alt", "Prob", "Mutation Prob.",
                "Exon Start", "Exon End", "Mutation Position", "Depth", "TumorDepth",
                "Peptide Normal", "Peptide Mutation", "DNA_Normal", "DNA_Mut",
                "TotalRNA", "TumorRNARatio", "TumorRNA", "nA", "nB", "Checker",
                "MutRatio", "MutRatio Min", "MutRatio Max")   

remove<-NULL
result<-NULL
type<-NULL
rna<-NULL
synonym<-NULL
for(file in files){
   #tumor and HLA
   fs<-list.files(paste(root, file, sep="/"))
   #HLA
   hit<-grep("HLAtype.txt", fs)
   f<-fs[hit]
   if(length(f)==0) {
      remove<-c(remove,file)
      next
   }
   hla<-scan(paste(root, file, f, sep="/"), "character", sep="\n")[3]
   #Synonym/Nonsynonym
   hit<-grep("summary.txt", fs)
   f<-fs[hit]
   if(length(f)==0) {
      remove<-c(remove,file)
      next
   }
   sy<-scan(paste(root, file, f, sep="/"), "character", sep="\n")
   num_sy<-as.numeric(strsplit(gsub("[ ]+", "\t",sy[match("exonic - synonymous", sy)+1]),"\t")[[1]][2])
   num_ny<-as.numeric(strsplit(gsub("[ ]+", "\t",sy[match("exonic - nonsynonymous", sy)+1]),"\t")[[1]][2])
   #Main
   hit<-grep("count.txt", fs)
   f<-fs[hit]
   if(length(f)==0) {
      remove<-c(remove,file)
      next
   }
   nn<-as.numeric(strsplit(scan(paste(root, file, f, sep="/"), "character", sep="\n")[4], "\t")[[1]][2])
   vh<-as.numeric(strsplit(scan(paste(root, file, f, sep="/"), "character", sep="\n")[5], "\t")[[1]][2])
   count<-sapply(scan(paste(root, file, f, sep="/"), "character", sep="\n")[6:43], function(x) strsplit(x, "\t")[[1]][2])
   count<-count[-c(13,26)]
   result<-cbind(result, c(as.numeric(count), 
                           as.numeric(count) / nn, 
			   as.numeric(count) / (vh * nn) * 6.0,
                           as.numeric(count) / num_sy, 
			   as.numeric(count) / (vh * num_sy) * 6.0))
   type<-c(type, hla)
   rna<-c(rna, as.numeric(strsplit(scan(paste(root, file, f, sep="/"), "character", sep="\n")[18], "\t")[[1]][2]))
   synonym<-c(synonym, num_ny / num_sy)
}

e<-c("min50HLAMer","min50HLASep","min500HLAMer","min500HLASep","min50wt500HLAMer","min50wt500HLASep",
     "rank0.5HLAMer","rank0.5HLASep","rank2HLAMer","rank2HLASep","rank0.5wt2HLAMer","rank0.5wt2HLASep")
e<-c(e, paste(e,"_fpkm_1",sep="_"),paste(e,"_fpkm_0",sep="_"))
rownames(result)<-c(e, paste(e,"/Nonsynonymous",sep="_"), 
		       paste(e,"/Nonsynonymous-HLA",sep="_"),
		       paste(e,"/Synonymous",sep="_"),
		       paste(e,"/Synonymous-HLA",sep="_"))

pdf(height = 24, width = 16, file = "SNVResult.pdf")
par(mai = c(0.85, 1.5, 0.68, 0.35))

synonym2<-synonym[which(synonym!=Inf)]
data<-sapply(unique(post), function(u) synonym2[which(!is.na(match(sapply(type[which(synonym!=Inf)], 
  				   function(y) strsplit(y,"-")[[1]][1]), pre[which(!is.na(match(post, u)))])))])
boxplot(data[order((as.numeric(lapply(data, median))))], horizontal=T, main="nonsy / synony", xlab="neoantigens / SNVs", las=1)

r<-1
for(i1 in 1:5){
 for(i2 in 1:(nrow(result)/(5*2))){
  data<-sapply(unique(post), function(u) result[r,which(!is.na(match(sapply(type, 
  			     function(y) strsplit(y,"-")[[1]][1]), pre[which(!is.na(match(post, u)))])))])
  boxplot(data[order((as.numeric(lapply(data, median))))], horizontal=T, main=rownames(result)[r], xlab="neoantigens / SNVs", las=1)
  r<-r+1
 }
 result_fpkm<-result[,which(as.numeric(rna)!=0)]
 if(i1==4 | i1==5) result_fpkm<-result[,which(as.numeric(rna)!=0 & synonym!=Inf)]
 for(i2 in 1:(nrow(result)/(5*2))){
  if(i1==4 | i1==5) {
    rem<-which(is.na(match(pre, sapply(type[which(as.numeric(rna)!=0 & synonym!=Inf)], function(x) strsplit(x, "-")[[1]][1]))))
    data<-sapply(unique(post[-rem]), function(u) result_fpkm[r, which(!is.na(match(sapply(type[which(as.numeric(rna)!=0 & synonym!=Inf)], 
  				   function(y) strsplit(y,"-")[[1]][1]),pre[which(!is.na(match(post, u)))])))])
  }else{
    rem<-which(is.na(match(pre, sapply(type[which(as.numeric(rna)!=0)], function(x) strsplit(x, "-")[[1]][1]))))
    data<-sapply(unique(post[-rem]), function(u) result_fpkm[r, which(!is.na(match(sapply(type[which(as.numeric(rna)!=0)], 
  				   function(y) strsplit(y,"-")[[1]][1]),pre[which(!is.na(match(post, u)))])))])
  }
  boxplot(data[order((as.numeric(lapply(data, median))))], horizontal=T, main=rownames(result)[r], xlab="neoantigens / SNVs", las=1)
  r<-r+1
 }
}

dev.off()
