#library(robustbase)
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
fpkm<-as.numeric(commandArgs(TRUE)[2])

index<-NULL
for(file in files){
   fs<-list.files(paste(dir, file, sep="/"))
   hit<-grep("HLAtype.txt", fs)
   f<-fs[hit]
   tumor<-sapply(scan(paste(dir, file, f, sep="/"), "character", sep="\n")[3], function(x) strsplit(x,"-")[[1]][1])
   hla<-paste(scan(paste(org, file, f, sep="/"), "character", sep="\n")[-c(1,2,3)], collapse=" ")
   index<-rbind(index, c(file, f, post[match(tumor,pre)], hla))
}
#RemoveHomo
index<-index[sapply(index[,4], function(x) length(unique(strsplit(x, " ")[[1]])))==6,]

remove<-NULL
#for(u in unique(index[,3])){
 u<-unique(index[,3])[as.numeric(commandArgs(TRUE)[1])]
 pdf(width=18,height=15,file=paste("StatSnvPersonMatrix/HLA.Boxplot.Class1.FULL.sp",u,fpkm,"pdf",sep="."))
 FULL<-NULL
 FULL_sep<-NULL
 HIT<-NULL
 tag<-NULL
 l<-which(!is.na(match(index[,3], u)))
 eachsubject<-NULL
 eachposition<-NULL
 for(i in l){
  #Original Data
  fs<-list.files(paste(org, index[i,1], sep="/"))
  hit<-grep("HLAtype.txt", fs)
  f<-fs[hit]
  if(length(f)==0) {
   remove<-c(remove,paste(org, index[i,1], sep="/"))
   next
  }
  hla<-scan(paste(org, index[i,1], f, sep="/"), "character", sep="\n")[2]
  hla<-paste(c(index[i,1], strsplit(hla," ")[[1]][-1]),collapse=" ")
  sp<-map[match(gsub(":","",strsplit(hla," ")[[1]][-1]),map[,1]),2]

  flag<-FALSE
  len<-1
  #IC50
  for(criteria in c("min_lower50.txt","min_upper500.txt")){
   fs<-list.files(paste(dir, index[i,1], sep="/"))
   hit<-grep("min_lower50.txt", fs)
   f<-fs[hit]
   if(length(f)==0) {
    remove<-c(remove,paste(dir, index[i,1], sep="/"))
    flag<-TRUE
    break
   }
   peptide<-t(sapply(scan(paste(dir, index[i,1], f, sep="/"), "character", sep="\n", skip=1), function(x) strsplit(x,"\t")[[1]]))
   peptide<-peptide[which(peptide[,21]!="NA"),]
   colnames(peptide)<-strsplit(scan(paste(dir, index[i,1], f, sep="/"), "character", sep="\n", nlines=1), "\t")[[1]]

   #Master
   if(fpkm!=0) {
     peptide<-peptide[!is.na(as.numeric(peptide[,28])) & as.numeric(peptide[,28]) > fpkm,]
   }
   if(nrow(peptide)==0) next
   for(ic50 in c(50, 200, 500)){
    tmp<-NULL
    tmp_sep<-NULL
    tmp2<-NULL
    for(h in 1:length(l)){
     hla_mer<-strsplit(index[l[h],4], " ")[[1]]
     count<-0
     count_sep<-0
     for(up in unique(peptide[,21])){
      hit<-which(!is.na(match(peptide[,21], up)))
      if(length(hit)==1){
       print(index[i,1])
       flag<-TRUE
       break
      }
      peptide_hit<-peptide[hit[sort(match(paste("HLA-", hla_mer, sep=""), peptide[hit,1]))],]
      if(min(as.numeric(peptide_hit[,5]))< ic50) count<-count+1.0/length(unique(peptide[,21]))
      count_sep<-count_sep + (length(which(as.numeric(peptide_hit[,5]) <ic50)) / nrow(peptide_hit)) / length(unique(peptide[,21]))
     }
     if(flag)break
     tmp<-c(tmp, count)
     tmp_sep<-c(tmp_sep, count_sep)
     tmp2<-c(tmp2, length(which(!is.na(match(hla_mer, strsplit(hla, " ")[[1]][-1])))))
    }
    if(flag)break
    if(length(FULL) < len){FULL<-c(FULL, list(NULL))}
    FULL[[len]]<-rbind(FULL[[len]], tmp)
    rownames(FULL[[len]])[nrow(FULL[[len]])]<-hla

    if(length(FULL_sep) < len){FULL_sep<-c(FULL_sep, list(NULL))}
    FULL_sep[[len]]<-rbind(FULL_sep[[len]], tmp_sep)
    rownames(FULL_sep[[len]])[nrow(FULL_sep[[len]])]<-hla

    if(length(HIT) < len){HIT<-c(HIT, list(NULL))}
    HIT[[len]]<-rbind(HIT[[len]], tmp2)
    rownames(HIT[[len]])[nrow(HIT[[len]])]<-hla

    tag<-c(tag, paste(u, criteria, ic50, sep="-"))
    len<-len+1
   }
  }
  if(flag) next

  #Rank
  for(criteria in c("rank_lower0.5.txt","rank_upper2.txt")){
   fs<-list.files(paste(dir, index[i,1], sep="/"))
   hit<-grep("min_lower50.txt", fs)
   f<-fs[hit]
   if(length(f)==0) {
    flag<-TRUE
    remove<-c(remove,paste(dir, index[i,1], sep="/"))
    break
   }
   peptide<-t(sapply(scan(paste(dir, index[i,1], f, sep="/"), "character", sep="\n", skip=1), function(x) strsplit(x,"\t")[[1]]))
   peptide<-peptide[which(peptide[,21]!="NA"),]
   colnames(peptide)<-strsplit(scan(paste(dir, index[i,1], f, sep="/"), "character", sep="\n", nlines=1), "\t")[[1]]
   if(fpkm!=0) {
     peptide<-peptide[!is.na(as.numeric(peptide[,28])) & as.numeric(peptide[,28]) > fpkm,]
   }

   for(ic50 in c(0.1, 0.5, 1, 2)){
    tmp<-NULL
    tmp_sep<-NULL
    tmp2<-NULL
    for(h in 1:length(l)){
     hla_mer<-strsplit(index[l[h],4], " ")[[1]]
     count<-0
     count_sep<-0
     for(up in unique(peptide[,21])){
      hit<-which(!is.na(match(peptide[,21], up)))
      if(length(hit)==1){
       print(index[i,1])
       flag<-TRUE
       break
      }
      peptide_hit<-peptide[hit[sort(match(paste("HLA-", hla_mer, sep=""), peptide[hit,1]))],]
      if(min(as.numeric(peptide_hit[,6]))< ic50) count<-count+1/length(unique(peptide[,21]))
      count_sep<-count_sep + (length(which(as.numeric(peptide_hit[,6]) <ic50)) / nrow(peptide_hit)) / length(unique(peptide[,21]))
     }
     if(flag)break
     tmp<-c(tmp, count)
     tmp_sep<-c(tmp_sep, count_sep)
     tmp2<-c(tmp2, length(which(!is.na(match(hla_mer, strsplit(hla, " ")[[1]][-1])))))
    }
    if(flag)break
    if(length(FULL) < len){FULL<-c(FULL, list(NULL))}
    FULL[[len]]<-rbind(FULL[[len]], tmp)
    rownames(FULL[[len]])[nrow(FULL[[len]])]<-hla

    if(length(FULL_sep) < len){FULL_sep<-c(FULL_sep, list(NULL))}
    FULL_sep[[len]]<-rbind(FULL_sep[[len]], tmp_sep)
    rownames(FULL_sep[[len]])[nrow(FULL_sep[[len]])]<-hla

    if(length(HIT) < len){HIT<-c(HIT, list(NULL))}
    HIT[[len]]<-rbind(HIT[[len]], tmp2)
    rownames(HIT[[len]])[nrow(HIT[[len]])]<-hla

    tag<-c(tag, paste(u, criteria, ic50, sep="-"))
    len<-len+1
   }
  }
  if(flag) next
 }

 print(length(FULL))
 for(j in 1:length(FULL)){
  print(dim(FULL[[j]]))
  if(ncol(FULL[[j]])==length(index[l,1])){
   FULL[[j]]<-FULL[[j]][,which(!is.na(match(index[l,1], sapply(rownames(FULL[[j]]), function(x) strsplit(x," ")[[1]][1]))))]
  }
  if(ncol(HIT[[j]])==length(index[l,1])){
   HIT[[j]]<-HIT[[j]][,which(!is.na(match(index[l,1], sapply(rownames(HIT[[j]]), function(x) strsplit(x," ")[[1]][1]))))]
  }
  val<-as.numeric(FULL[[j]])
  x<-as.numeric(HIT[[j]])
  boxplot(sapply(0:6, function(y) val[grep(y, x)]), names=0:6, xlab="Num. Hit", 
  		      ylab="Ave.Num.NeoAntigen",main=tag[j])

  #Mizuno
  Ftmp<-t(apply(FULL[[j]], 1, scale))
  Ftmp[which(is.na(apply(Ftmp, 1, sum))),]<-FULL[[j]][which(is.na(apply(Ftmp, 1, sum))),]
  less<-length(which(sapply(1:ncol(Ftmp), function(x) length(which(Ftmp[x,x] > Ftmp[,x])) < length(which(Ftmp[x,x] < Ftmp[,x])))))
  more<-length(which(sapply(1:ncol(Ftmp), function(x) length(which(Ftmp[x,x] > Ftmp[,x])) > length(which(Ftmp[x,x] < Ftmp[,x])))))
  boxplot(Ftmp, xlim=c(0.5,ncol(Ftmp)+0.5), ylim=c(0,max(Ftmp)*1.1), main="Fix HLAs",
   xlab=paste("Subject", less, more, round(binom.test(c(less,more),0.5)$p.val, 5)), ylab="Ratio")
  par(new=TRUE)
  plot(1:ncol(Ftmp), sapply(1:ncol(Ftmp), function(x) Ftmp[x,x]),pch=4,
   xlim=c(0.5,ncol(Ftmp)+0.5),ylim=c(0,max(Ftmp)*1.1),xlab="", ylab="")

  Ftmp<-FULL[[j]]
  less<-length(which(sapply(1:nrow(Ftmp), function(x) length(which(Ftmp[x,x] > Ftmp[x,])) < length(which(Ftmp[x,x] < Ftmp[x,])))))
  more<-length(which(sapply(1:nrow(Ftmp), function(x) length(which(Ftmp[x,x] > Ftmp[x,])) > length(which(Ftmp[x,x] < Ftmp[x,])))))
  boxplot(t(Ftmp), xlim=c(0.5,nrow(Ftmp)+0.5), ylim=c(0,max(Ftmp)*1.1), main="Fix Variants",
   xlab=paste("Subject", less, more, round(binom.test(c(less,more),0.5)$p.val, 5)), ylab="Ratio")
  par(new=TRUE)
  plot(1:nrow(Ftmp), sapply(1:nrow(Ftmp), function(x) Ftmp[x,x]),pch=4,
   xlim=c(0.5,nrow(Ftmp)+0.5),ylim=c(0,max(Ftmp)*1.1),xlab="", ylab="")

  eachposition<-rbind(eachposition, cbind(tag[j], sapply(rownames(FULL[[j]]), function(x) strsplit(x," ")[[1]][1]),
   sapply(1:nrow(FULL[[j]]), function(x) length(which(FULL[[j]][x,x] > FULL[[j]][x,]))) / ncol(FULL[[j]])))

  dist=NULL
  y_range<-c(0,1)
  count<-1
  ylim_<-c(max(0, min(val) * 0.9), min(1, max(val)* 1.1))
  plot(0:6, sapply(0:6, function(z) median(FULL[[j]][1,grep(z, HIT[[j]][1,])])), 
   type="b", xlim=c(0,6), ylim=ylim_, col=count, pch=count)
  dist<-c(dist, lm(FULL[[j]][1, ] ~ HIT[[j]][1, ])$coefficients[2])
  abline(lm(FULL[[j]][1, ] ~ HIT[[j]][1, ]), col=count, pch=count, lty=3)
  for(k in 1:nrow(FULL[[j]])){
   count<-count+1
   par(new=TRUE)
   plot(0:6, sapply(0:6, function(z) median(FULL[[j]][k, grep(z, HIT[[j]][k,])])), 
    type="b", xlim=c(0, 6), ylim=ylim_, xlab="", ylab="", xaxt="n", yaxt="n", col=count, pch=count)
   dist<-c(dist, lm(FULL[[j]][k, ] ~ HIT[[j]][k, ])$coefficients[2])
   abline(lm(FULL[[j]][k, ] ~ HIT[[j]][k, ]), col=count, pch=count, lty=3)
   eachsubject<-rbind(eachsubject, c(tag[j], strsplit(rownames(FULL[[j]])[k]," ")[[1]][1], lm(FULL[[j]][k, ] ~ HIT[[j]][k, ])$coefficients[2]))
  }
  hist(dist, breaks="Scott", main=paste("pval=", t.test(dist)$p.value))
 }

 print(length(FULL_sep))
 for(j in 1:length(FULL_sep)){
  print(dim(FULL_sep[[j]]))
  if(ncol(FULL_sep[[j]])==length(index[l,1])){
   FULL_sep[[j]]<-FULL_sep[[j]][,which(!is.na(match(index[l,1], sapply(rownames(FULL_sep[[j]]), function(x) strsplit(x," ")[[1]][1]))))]
  }
  val<-as.numeric(FULL_sep[[j]])
  x<-as.numeric(HIT[[j]])
  boxplot(sapply(0:6, function(y) val[grep(y, x)]), names=0:6, xlab="Num. Hit", 
  		      ylab="Ave.Num.NeoAntigen",main=paste("Sep", tag[j]))

  #Mizuno
  Ftmp<-t(apply(FULL_sep[[j]], 1, scale))
  Ftmp[which(is.na(apply(Ftmp, 1, sum))),]<-FULL_sep[[j]][which(is.na(apply(Ftmp, 1, sum))),]
  less<-length(which(sapply(1:ncol(Ftmp), function(x) length(which(Ftmp[x,x] > Ftmp[,x])) < length(which(Ftmp[x,x] < Ftmp[,x])))))
  more<-length(which(sapply(1:ncol(Ftmp), function(x) length(which(Ftmp[x,x] > Ftmp[,x])) > length(which(Ftmp[x,x] < Ftmp[,x])))))
  boxplot(Ftmp, xlim=c(0.5,ncol(Ftmp)+0.5), ylim=c(0,max(Ftmp)*1.1), main="Fix HLAs",
   xlab=paste("Subject", less, more, round(binom.test(c(less,more),0.5)$p.val, 5)), ylab="Ratio")
  par(new=TRUE)
  plot(1:ncol(Ftmp), sapply(1:ncol(Ftmp), function(x) Ftmp[x,x]),pch=4,
   xlim=c(0.5,ncol(Ftmp)+0.5),ylim=c(0,max(Ftmp)*1.1),xlab="", ylab="")

  Ftmp<-FULL_sep[[j]]
  less<-length(which(sapply(1:nrow(Ftmp), function(x) length(which(Ftmp[x,x] > Ftmp[x,])) < length(which(Ftmp[x,x] < Ftmp[x,])))))
  more<-length(which(sapply(1:nrow(Ftmp), function(x) length(which(Ftmp[x,x] > Ftmp[x,])) > length(which(Ftmp[x,x] < Ftmp[x,])))))
  boxplot(t(Ftmp), xlim=c(0.5,nrow(Ftmp)+0.5), ylim=c(0,max(Ftmp)*1.1), main="Fix Variants",
   xlab=paste("Subject", less, more, round(binom.test(c(less,more),0.5)$p.val, 5)), ylab="Ratio")
  par(new=TRUE)
  plot(1:nrow(Ftmp), sapply(1:nrow(Ftmp), function(x) Ftmp[x,x]),pch=4,
   xlim=c(0.5,nrow(Ftmp)+0.5),ylim=c(0,max(Ftmp)*1.1),xlab="", ylab="")

  eachposition<-rbind(eachposition, cbind(paste("Sep", tag[j]), sapply(rownames(FULL_sep[[j]]), function(x) strsplit(x," ")[[1]][1]),
   sapply(1:nrow(FULL_sep[[j]]), function(x) length(which(FULL_sep[[j]][x,x] > FULL_sep[[j]][x,]))) / ncol(FULL_sep[[j]])))

  dist=NULL
  y_range<-c(0,1)
  count<-1
  ylim_<-c(max(0, min(val) * 0.9), min(1, max(val)* 1.1))
  plot(0:6, sapply(0:6, function(z) median(FULL_sep[[j]][1,grep(z, HIT[[j]][1,])])), 
   type="b", xlim=c(0,6), ylim=ylim_, col=count, pch=count)
  dist<-c(dist, lm(formula = FULL_sep[[j]][1, ] ~ HIT[[j]][1, ])$coefficients[2])
  abline(lm(formula = FULL_sep[[j]][1, ] ~ HIT[[j]][1, ]), col=count, pch=count, lty=3)
  for(k in 1:nrow(FULL_sep[[j]])){
   count<-count+1
   par(new=TRUE)
   plot(0:6, sapply(0:6, function(z) median(FULL_sep[[j]][k, grep(z, HIT[[j]][k,])])), 
    type="b", xlim=c(0, 6), ylim=ylim_, xlab="", ylab="", xaxt="n", yaxt="n", col=count, pch=count)
   dist<-c(dist, lm(FULL_sep[[j]][k, ] ~ HIT[[j]][k, ])$coefficients[2])
   abline(lm(FULL_sep[[j]][k, ] ~ HIT[[j]][k, ]), col=count, pch=count, lty=3)
   eachsubject<-rbind(eachsubject, c(paste("Sep", tag[j]), strsplit(rownames(FULL[[j]])[k]," ")[[1]][1], lm(FULL[[j]][k, ] ~ HIT[[j]][k, ])$coefficients[2]))
  }
  hist(dist, breaks="Scott", main=paste("pval=", t.test(dist)$p.value))
 }
#}

write.table(eachsubject,paste("StatSnvPersonMatrix/HLA.Boxplot.Class1.FULL.sp.lm",u,fpkm,"txt",sep="."),
        row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")
write.table(eachposition,paste("StatSnvPersonMatrix/HLA.Boxplot.Class1.FULL.sp.position",u,fpkm,"txt",sep="."),
        row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")
print("Remove")
print(remove)
dev.off()

