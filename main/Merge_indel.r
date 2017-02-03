dir1<-commandArgs(TRUE)[1]
fols<-list.files(commandArgs(TRUE)[1])
start<-as.numeric(commandArgs(TRUE)[2])
end<-as.numeric(commandArgs(TRUE)[3])
for(fol in fols[start:end]){
   print(fol)
   dir2<-paste(dir1,fol, sep="/")
   files<-list.files(paste(dir1,fol, sep="/"))
   if(length(grep("count",files))==1){
      next
   }
   #Get Peptide Info
   files_part<-files[is.na(sapply(files, function(x) strsplit(x, "peptide.indel.txt")[[1]][2]))]
   file<-files_part[grep("peptide.indel.txt", files_part)]
   file<-file[order(sapply(file, nchar))[1]]
   if(is.na(file)){
      print("No File:!!")
      next
   }
   prefix<-file
   if(length(scan(paste(dir2, file, sep="/"), "character", sep="\n"))==0){
      next
   }
   info<-t(sapply(scan(paste(dir2, file, sep="/"), "character", sep="\n"), 
				       function(x) strsplit(x, "\t")[[1]]))
   cinfo<-c("", "Gene ID", "Chr", "NM_ID", "Change", "ref", "alt", "Prob", "Mutation Prob.", 
   		"Exon Start", "Exon End", "Mutation Position", "Depth", "TumorDepth", 
		"Peptide Normal", "Peptide Mutation", "DNA_Normal", "DNA_Mut", 
		"TotalRNA", "TumorRNARatio", "TumorRNA", "nA", "nB", "Checker", 
		"MutRatio", "MutRatio Min", "MutRatio Max")
   info<-info[,1:length(cinfo)]
   if(is.null(ncol(info))){info<-t(as.matrix(info))}
   colnames(info)<-cinfo
   
   #Remove RNAseq Info
   rownames(info)<-NULL
   info<-info[,-match(c("DNA_Normal", "DNA_Mut"), colnames(info))]
   if(is.null(ncol(info))){info<-t(as.matrix(info))}

   #Get Peptide
   FASTA<-files[grep("\\.fasta\\.HLACLASS", files)]
   FASTA<-FASTA[grep("\\.peptide.indel.txt", FASTA)]
   if(length(FASTA)==0) {
      if(length(scan(paste(dir2,files[grep("HLAtype.txt",files)],sep="/"), "character", sep="\n"))!=3){
         print("Error:")
      }
      next
   }

   e<-NULL
   #FASTA_list<-c(as.list(FASTA),list(FASTA))
  #for(FASTA in FASTA_list){
   full_peptide<-NULL
   min_peptide<-NULL
   rank_peptide<-NULL
   HLA_count<-0
   for(f in FASTA){
      print(paste(dir2, f, sep="/"))
      test1<-scan(paste(dir2, f, sep="/"), "character", sep="\n", skip=1)
      test1<-gsub(" <= WB", "", test1)
      ss1<-grep(" Pos ", test1)+2
      ee1<-grep("Protein", test1)-2
      num1<-sapply(gsub("[ ]+","\t",test1[ss1]), function(x) strsplit(x, "\t")[[1]][12])

      if(length(grep("cannot be found in hla_pseudo list", test1))>0) next
      HLA_count <- HLA_count + 1
      for(h1 in 1:length(num1)){
         d1<-t(sapply(gsub("[ ]+", "\t", test1[ss1[h1]:ee1[h1]]), function(x) strsplit(x, "\t")[[1]][c(3,2,12,11,4,14,15)]))
	 if(nrow(d1)==0) {
	    r_can<-match(num1[h1], info[,2])
	    if(is.na(r_can)){r_can<-grep(num1[h1], info[,2])}
	    remove<-c(remove, r_can)
	    next
	 }
	 full_peptide<-rbind(full_peptide, d1)
      }
   }
   
   if(is.null(full_peptide)) next
   if(nrow(full_peptide)==0) next
   tag<-c("HLA", "Pos", "Gene", "Core", "MutatedPeptide", "Mut_IC50", "Mut_Rank", colnames(info))
   remove<-which(sapply(full_peptide[,5], function(x) length(grep(x, min_peptide[1,22])))==1)
   if(length(remove)>0) full_peptide<-full_peptide[-remove, ]

   #Bind Peptide and info
   full_peptide<-cbind(full_peptide, info[match(full_peptide[,3], info[,2]),])
   colnames(full_peptide)<-tag
   write.table(full_peptide, paste(dir2, "/", prefix, ".ALL.txt", sep=""), row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")

   tmp<-NULL
   if(nrow(full_peptide) == 1){
    full_peptide<-t(full_peptide)
   }
   for(HLA in c(list(unique(full_peptide[,1])),unique(full_peptide[,1]))){
    min_peptide_h<-full_peptide[which(!is.na(match(full_peptide[,1], HLA))),]
    for(gene in unique(full_peptide[,3])){
     num_list<-NULL
     for(th in c(25, 50, 100, 200, 500)){
      num<-0
      min_peptide<-min_peptide_h[which(!is.na(match(min_peptide_h[,3], gene))),]
      for(uqh in unique(min_peptide[,4])){
       index<-which(!is.na(match(min_peptide[,4], uqh)))
       if(min(as.numeric(min_peptide[index,6])) < th) num<-num+1
      }
      num_list<-c(num_list, num)
     }
     for(th in c(0.5, 1, 2, 5)){
      num<-0
      min_peptide<-min_peptide_h[which(!is.na(match(min_peptide_h[,3], gene))),]
      for(uqh in unique(min_peptide[,4])){
       index<-which(!is.na(match(min_peptide[,4], uqh)))
       if(min(as.numeric(min_peptide[index,7])) < th) num<-num+1
      }
      num_list<-c(num_list, num)
     }
     HLA_name<-"full"
     if(length(HLA)==1)HLA_name<-HLA
     tmp<-rbind(tmp, c(HLA_name, gene, num_list, nchar(min_peptide[1,23]) - nchar(min_peptide[1,22]), min_peptide[1,24], min_peptide[1,25], min_peptide[1,26],
    	       min_peptide[1,30],min_peptide[1,31],min_peptide[1,32]))
    }
   }
   colnames(tmp)<-c("HLA","Gene","25","50","100","200","500","0.5","1","2","5","Length","TotalRNA","Depth","TumorRNA","CCFP","CCFP_min","CCFP_max")
   write.table(tmp, paste(dir2, "/", prefix, ".con.txt", sep=""), row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")
}

