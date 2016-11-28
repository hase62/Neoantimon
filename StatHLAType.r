data1<-sapply(scan("release_aug2015.v1.tsv.cut.txt", "character", sep="\n", skip=1), function(x) strsplit(x, "\t")[[1]][1])
data2<-scan("r1.txt", "character", sep="\n",skip=1)
data3<-scan("AugPilot_HLAGenotype.tsv", "character", sep="\n")
data4<-sapply(data3, function(x) strsplit(x, "\t")[[1]][1])
data5<-data2[match(data4, data1)]

write.table("", "HLATypesEachGroup.txt", row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t",append=FALSE)
for(u in unique(data5)){
   group<-which(!is.na(match(data5,u)))
   freq<-paste(sapply(data3[which(!is.na(match(data5,u)))],function(x) strsplit(x,"\t")[[1]][-1]))   
   freq<-strsplit(freq, "\"")
   hla<-sort(unique(unlist(lapply(freq, function(x) x[nchar(x)>5]))))
   hla<-hla[nchar(hla)==11]
   write.table(t(c(u, hla)), "HLATypesEachGroup.txt", row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t",append=TRUE)
}

doner<-sapply(scan("release_aug2015.v1.tsv.cut.txt", "character", sep="\n", skip=1), function(x) strsplit(x, "\t")[[1]][2])
write.table("", "DonerNameEachGroup.txt", row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t",append=FALSE)
for(u in unique(data5)){
   write.table(t(c(u, doner[which(!is.na(match(data5, u)))])), 
      "DonerNameEachGroup.txt", row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t",append=TRUE)
}

