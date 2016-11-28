library("meta")
library("rmeta")

org<-"Output.Snv.1"
dir<-"Output.HLA.1"
files<-list.files(dir)

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

map<-t(sapply(scan("SuperClass.txt","character",sep="\n"), function(x) strsplit(x,"\t")[[1]]))
#If you do not use SuperClass, 
#map[,1]<-"XXX"

index<-NULL
for(file in files){
   fs<-list.files(paste(dir, file, sep="/"))
   hit<-grep("HLAtype.txt", fs)
   f<-fs[hit]
   tumor<-sapply(scan(paste(dir, file, f, sep="/"), "character", sep="\n")[3], function(x) strsplit(x,"-")[[1]][1])
   hla<-scan(paste(dir, file, f, sep="/"), "character", sep="\n")[-c(1,2,3)]
   if(length(hla)==3) {
      next
   }
   index<-rbind(index, c(file, f, post[match(tumor,pre)], length(hla)))
}

remove<-NULL
Sum<-NULL
pdf(width=18, height=15,file="MetaPlot.Class1.sp.pdf")
for(u in unique(index[,3])){
 FULL<-NULL
 l<-which(!is.na(match(index[,3], u)))
 hla_list<-NULL
 #Make Matrixes for Each Tumor
 for(i in l){
  #Original Data
  fs<-list.files(paste(org, index[i,1], sep="/"))
  hit<-grep("HLAtype.txt", fs)
  f<-fs[hit]
  if(length(f)==0) {
   remove<-c(remove,file)
   hla_list<-hla_list[-length(hla_list)]
   next
  }
  hla<-scan(paste(org, index[i,1], f, sep="/"), "character", sep="\n")
  if(length(hla)<(3+6)) {
   remove<-c(remove,file)
   hla_list<-hla_list[-length(hla_list)]
   next
  }
  hla<-hla[2]
  sp<-map[match(gsub(":","",strsplit(hla," ")[[1]][-1]),map[,1]),2]
  hla_list<-c(hla_list, 
              paste(ifelse(is.na(sp)|sp=="Unclassified",gsub(":","",strsplit(hla," ")[[1]][-1]),sp),collapse=","))

  #HLA Shuffle
  fs<-list.files(paste(dir, index[i,1], sep="/"))
  hit<-grep("count.txt", fs)
  f<-fs[hit]
  if(length(f)==0) {
   remove<-c(remove,file)
   hla_list<-hla_list[-length(hla_list)]
   next
  }
  data<-t(sapply(scan(paste(dir, index[i,1], f, sep="/"), "character", sep="\n"),
        function(x) strsplit(x,"\t")[[1]]))
  rownames(data)<-NULL
  colnames(data)<-NULL
  data<-data[,-ncol(data)]
  rownames(data)<-data[,1]
  data<-data[,-1]
  skip<-grep("Length", rownames(data))
  num<-6:nrow(data)
  rdata<-NULL
  ny<-as.numeric(data[4,1])
  for(j in num){
   if(length(grep("Length", rownames(data)[j]))>0) {
    ny<-as.numeric(data[j,1])
    next
   }
   if(ny==0) ny<-0.1
   rdata<-rbind(rdata, as.numeric(data[j,])/ny)
   rownames(rdata)[nrow(rdata)]<-rownames(data)[j]
  }
  colnames(rdata)<-data[2,]

  if(!is.null(FULL)) {if(ncol(rdata)!=ncol(FULL[[1]])){
   remove<-c(remove,file)
   hla_list<-hla_list[-length(hla_list)]
   next
  }}
  FULL<-c(FULL,list(rdata))
 }

 #MetaAnalysis
 if(length(FULL) < 3) next
 for(i in 1:nrow(FULL[[1]])){
  meta<-NULL
  for(j in 1:length(FULL)){
   meta<-rbind(meta, FULL[[j]][i,])
  }
  nonzero<-which(!is.na(apply(meta, 1, sum)))
  meta<-meta[nonzero, ]
  sp<-map[match(gsub(":","",sapply(colnames(meta), function(x) strsplit(x,"-")[[1]][2])),map[,1]),2]
  hla2<-colnames(meta)<-ifelse(is.na(sp)|sp=="Unclassified",gsub(":","",sapply(colnames(meta), 
  				   function(x) strsplit(x,"-")[[1]][2])),sp)
  Bin<-c(0,0)
  tres<-NULL
  for(h in 1:length(unique(hla2))){
   tmp<-sapply(hla_list[nonzero], function(x) strsplit(x,",")[[1]])
   #Person who has unique(hla2)[h]
   if(is.list(tmp))   on<-which(unlist(lapply(tmp, function(x) !is.na(match(unique(hla2)[h], x)))))
   if(is.matrix(tmp)) on<-which(unlist(apply(tmp,2,function(x) !is.na(match(unique(hla2)[h], x)))))
   if(length(on) < 3) next
   mean_a<-meta[on,h]
   mean_b<-meta[-on,h]
   if(mean(mean_a) < mean(mean_b)) {Bin[1]<-Bin[1]+1
   }else{Bin[2]<-Bin[2]+1}
   tres<-rbind(tres, c(t.test(mean_a,mean_b)$p.value,length(mean_a),mean(mean_a),sqrt(var(mean_a)),length(mean_b),mean(mean_b),sqrt(var(mean_b))))
   rownames(tres)[nrow(tres)]<-unique(hla2)[h]
  }
  tres<-tres[!is.na(tres[,1]),]
  if(is.null(nrow(tres))) next
  if(nrow(tres) < 3) next
  tres[,4]<-ifelse(tres[,4]==0, tres[,7] * 0.01, tres[,4])
  forest(metacont(n.e=tres[,2], mean.e=tres[,3], sd.e=tres[,4], n.c=tres[,5], mean.c=tres[,6], sd.c=tres[,7], 
   sm="SMD",studlab=rownames(tres)),hetlab=paste(u,rownames(FULL[[1]])[i]))
  Sum<-rbind(Sum, Bin)
  rownames(Sum)[nrow(Sum)]<-paste(u,rownames(FULL[[1]])[i])
 }
}
dev.off()
write.table(remove, "RmList.txt", row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")
write.table(Sum, "StatSnvMatrix.class1.sp.txt", row.names=TRUE, col.names=FALSE, quote=FALSE, sep="\t")
