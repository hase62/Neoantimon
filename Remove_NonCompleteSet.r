output<-commandArgs(TRUE)[1]
files<-list.files(output)
rm_list<-NULL
for(f in files){
   fs<-list.files(paste(output, f, sep="/"))

   if(length(fs)==0){
     rm_list<-c(rm_list, paste(output,f,sep="/"))
     next
   }

   #If no estiamated HLA
   HLAs<-strsplit(scan(paste(output,f,fs[grep("HLAtype",fs)],sep="/"),"character",sep="\n")[2], " ")[[1]][-1]
   uHLAs<-unique(sapply(HLAs, function(x) strsplit(x,"\\*")[[1]][1]))
   if(length(HLAs)==0|(length(grep("DP",uHLAs))==1 & length(grep("DQ",uHLAs))==1 & length(grep("DR",uHLAs))==0)) {
      print("HLA not determined.")
      print(f)
      next
   }
   if(length(fs) == 4|length(fs) == 8) {
     if(length(grep("\\.extracted\\.txt", fs))>0){
      fs<-fs[grep("\\.extracted\\.txt", fs)]
      if(length(grep("\\.extracted\\.txt\\.", fs))>0) {
       fs<-fs[-grep("\\.extracted\\.txt\\.", fs)]
      }
      if(length(scan(paste(output,f,fs[1],sep="/"),"character",sep="\n"))==1) {
       next
      }
     }
     rm_list<-c(rm_list, paste(output,f,sep="/"))
     next
   }
   if(length(fs) < 20){
     if(length(grep("\\.extracted\\.txt", fs))>0){
      fss<-fs[grep("\\.extracted\\.txt", fs)]
      if(length(grep("\\.extracted\\.txt\\.", fss))>0) {
       fss<-fss[-grep("\\.extracted\\.txt\\.", fss)]
      }
      if(length(scan(paste(output,f,fss[1],sep="/"),"character",sep="\n"))==1) {
       next
      }
     }
     if(length(grep("HLACLASS",fs))==0){
      rm_list<-c(rm_list, paste(output,f,sep="/"))
      next
     }
   }
   if(length(grep("HLAtype", fs)) < 1) {
      rm_list<-c(rm_list, paste(output,f,sep="/"))
      next
   }
   if(length(grep("indel",fs))>0){
      fs_temp<-fs[grep("output.tsv.peptide.indel.txt",fs)]
   }else{
      fs_temp<-fs[grep("extracted.txt.peptide.txt",fs)]
   }

   nums<-as.numeric(unique(gsub("\\.","",sapply(fs[grep("fasta.HLACLASS",fs)], 
    function(x) strsplit(x,"HLACLASS1|HLACLASS2|peptide|norm")[[1]][3]))))
   if(length(nums)==0) {
      rm_list<-c(rm_list, paste(output,f,sep="/"))
      next
   }
   for(num in 1:max(nums)){
      if(length(grep("indel",fs))>0){
       if(length(grep("HLACLASS1",fs))>0){
        hit1<-grep(paste("\\.",num,"\\.HLACLASS1.",num,".peptide.indel.sh", sep=""), fs)
       }else{
        hit1<-grep(paste("\\.",num,"\\.HLACLASS2.",num,".peptide.indel.sh", sep=""), fs)
       }
      }else{
       hit1<-c(grep(paste("\\.",num,"\\.peptide.sh", sep=""), fs), grep(paste("\\.",num,"\\.normpeptide.sh", sep=""), fs))
      }
      if(length(hit1)%%2!=0 | length(hit1)==0) {
        rm_list<-c(rm_list, paste(output,f,sep="/"))
        break
      }
      if(length(grep("indel",fs))>0){
       if(length(grep("HLACLASS1",fs))>0){
        hit2<-grep(paste("HLACLASS1.",num,".peptide.indel",sep=""), fs)
       }else{
        hit2<-grep(paste("HLACLASS2.",num,".peptide.indel",sep=""), fs)
       }
      }else{
       if(length(grep("HLACLASS1",fs))>0){
        hit2<-c(grep(paste("HLACLASS1.",num,".peptide",sep=""), fs),grep(paste("HLACLASS1.",num,".normpeptide",sep=""), fs))
       }else{
        hit2<-c(grep(paste("HLACLASS2.",num,".peptide",sep=""), fs),grep(paste("HLACLASS2.",num,".normpeptide",sep=""), fs))
       }
      }
      if(length(hit2)%%2!=0 | length(hit1) == length(hit2)){
       rm_list<-c(rm_list, paste(output,f,sep="/"))
       break
      }
      if(length(hit1)==4){
       rm_list<-c(rm_list, paste(output,f,sep="/"))
       break
      }
   }
}
rm_list<-unique(rm_list)

write.table(rm_list, paste("RmList",output,"txt",sep="."), col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\t")

