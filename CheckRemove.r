dir<-commandArgs(TRUE)[1]
files<-list.files(dir)

rm<-NULL
for(file in files){
 fs<-list.files(paste(dir, file, sep="/"))
 if(length(grep("HLACLASS2.normpeptide.txt", fs))==0|length(grep("HLACLASS2.peptide.txt", fs))==0){
  rm<-c(rm,file)
  next
 }
 d1<-scan(paste(dir, file, fs[grep("HLACLASS2.normpeptide.txt", fs)[1]], sep="/"), "character", sep="\n")
 d2<-scan(paste(dir, file, fs[grep("HLACLASS2.peptide.txt", fs)[1]], sep="/"), "character", sep="\n")
 if(length(d1) < length(d2) * 0.95 | length(d1) > length(d2) * 1.05) rm<-c(rm, file)
 #if(length(d1) != length(d2)) rm<-c(rm, file)

}

write.table(paste(dir, rm,sep="/"), "RmList.txt", row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")
