dir1<-commandArgs(TRUE)[1]
fols<-list.files(commandArgs(TRUE)[1])
start<-as.numeric(commandArgs(TRUE)[2])
end<-as.numeric(commandArgs(TRUE)[3])

RemoveRNAList<-scan("RNABlackList.txt","character",sep="\n")
for(fol in fols[start:end]){
   print(fol)
   dir2<-paste(dir1,fol, sep="/")
   files<-list.files(paste(dir1,fol, sep="/"))
   if(length(grep("count",files))==1){
      next
   }
   #Get Peptide Info
   files_part<-files[is.na(sapply(files, 
                function(x) strsplit(x, "extracted.txt.peptide.txt")[[1]][2]))]
   file<-files_part[grep("extracted.txt.peptide.txt", files_part)]
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
   info[,match("TumorRNA",colnames(info))]<-as.numeric(info[,match("TotalRNA",colnames(info))]) * 
     as.numeric(info[,match("TumorDepth",colnames(info))]) / 
       as.numeric(info[,match("Depth",colnames(info))])   
   info[which(!is.na(match(sapply(info[,match("Gene ID",colnames(info))], 
     function(x) strsplit(x, "_")[[1]][2]), RemoveRNAList))),match("TumorRNA",colnames(info))]<-NA

   #Remove if PseudoGene
   disable<-NULL
   if(length(grep("Pse",commandArgs(TRUE)[1]))>0){
    ext<-scan(paste(dir2,gsub(".peptide.txt","",file),sep="/"),"character",sep="\n")
    disable<-sapply(ext[unique(c(grep("Disablements=0",ext),grep("AfterStop=0",ext)))],function(x) strsplit(x,"\t")[[1]][2])
   }

   #Remove RNAseq Info
   rownames(info)<-NULL
   info<-info[,-match(c("DNA_Normal", "DNA_Mut"), colnames(info))]
   if(is.null(ncol(info))){info<-t(as.matrix(info))}

   #Get Peptide
   FASTA<-files[grep("\\.fasta\\.HLACLASS", files)]
   FASTA<-FASTA[grep("\\.peptide.txt", FASTA)]
   if(length(FASTA)==0) {
      if(length(scan(paste(dir2,files[grep("HLAtype.txt",files)],sep="/"), "character", sep="\n"))!=3){
         print("Error:")
      }
      next
   }

   #Include Stop Codon
   remove<-which(sapply(info[,c(16)], function(x) length(grep("X", rev(strsplit(x, "")[[1]])[-1]))>0))
   #Remove XY
   remove<-c(remove, grep("X|Y",info[,3]))
   #Remove from Info
   if(length(remove) > 0) info<-info[-remove,]
   if(is.null(ncol(info))){info<-t(as.matrix(info))}
   if(nrow(info)==0) next

   e<-NULL
   FASTA_list<-c(as.list(FASTA),list(FASTA))
  for(FASTA in FASTA_list){
   full_peptide<-NULL
   min_peptide<-NULL
   min_peptide_50<-NULL
   rank_peptide<-NULL
   rank_peptide_50<-NULL
   HLA_count<-0
   #allele,start,end,length,peptide,ic50,Rank,Peptide_Normal_Sep,norm_ic_50,norm_Rank
   for(f in FASTA){
      print(paste(dir2, f, sep="/"))
      test1<-scan(paste(dir2, f, sep="/"), "character", sep="\n", skip=1)
      test1<-gsub(" <= WB", "", test1)
      ss1<-grep(" Pos ", test1)+2
      ee1<-grep("Protein", test1)-2
      num1<-sapply(gsub("[ ]+","\t",test1[ss1]), function(x) strsplit(x, "\t")[[1]][12])

      test2<-scan(paste(dir2, sub("peptide\\.txt", "normpeptide\\.txt", f), sep="/"),"character", sep="\n",skip=1)
      test2<-gsub(" <= WB", "", test2)
      ss2<-grep(" Pos ", test2)+2
      ee2<-grep("Protein", test2)-2
      num2<-sapply(gsub("[ ]+","\t",test2[ss2]), function(x) strsplit(x, "\t")[[1]][12])

      if(length(grep("cannot be found in hla_pseudo list", test1))>0) next
      HLA_count <- HLA_count + 1
      for(h1 in 1:length(num1)){
         #Skip if not match
	 if(is.na(grep(num1[h1], info[,2])[1]))next

         d4<-NULL
         hit<-match(num1[h1], num2)
         d1<-t(sapply(gsub("[ ]+", "\t", test1[ss1[h1]:ee1[h1]]), function(x) strsplit(x, "\t")[[1]][c(2,3,4,12,14,15)]))
      	 d2<-t(sapply(gsub("[ ]+", "\t", test2[ss2[hit]:ee2[hit]]), function(x) strsplit(x, "\t")[[1]][c(2,3,4,12,14,15)]))
	 l1<-sapply(d1[,3], nchar)
	 l2<-sapply(d2[,3], nchar)
	 for(r1 in unique(l1)){
	   hit1<-which(l1==r1)
	   hit2<-which(l2==r1)
	   if(length(hit1)==0) next
	   if(length(hit1)==1){
	    d3<-t(c(d1[hit1,c(2,1,4,3,5,6)], d2[hit2[match(d1[hit1,1], d2[hit2,1])],c(3,5,6)]))
	   }else{
	    d3<-cbind(d1[hit1,c(2,1,4,3,5,6)], d2[hit2[match(d1[hit1,1], d2[hit2,1])],c(3,5,6)])
	   }
	   d3<-d3[d3[,4]!=d3[,7],]
	   d4<-rbind(d4, d3)
	 }
	 if(nrow(d4)==0) {
	    print("Warning!! d4 is zero!!")
	    next
	 }
	 full_peptide<-rbind(full_peptide, d4)
	 min_peptide<-rbind(min_peptide, d4[order(as.numeric(d4[,5]))[1],])
	 rank_peptide<-rbind(rank_peptide, d4[order(as.numeric(d4[,6]))[1],])

	 d5<-d4[as.numeric(d4[,8]) > 500,]
	 if(is.vector(d5)){d5<-t(d5)}
	 if(nrow(d5)>0){
	  min_peptide_50<-rbind(min_peptide_50, d5[order(as.numeric(d5[,5]))[1],])
         }
	 d6<-d4[as.numeric(d4[,9]) > 2,]
	 if(is.vector(d6)){d6<-t(d6)}
	 if(nrow(d6)>0){
	  rank_peptide_50<-rbind(rank_peptide_50, d6[order(as.numeric(d6[,6]))[1],])
         }
      }
   }
   
   if(is.null(full_peptide)) next
   if(nrow(full_peptide)==0) next
   tag<-c("HLA", "Pos", "Gene", "MutatedPeptide", "Mut_IC50", "Mut_Rank", 
          "Norm_Peptide", "Norm_IC50", "Norm_Rank", colnames(info))

   #Bind Full Peptide and info
   if(nrow(full_peptide)==1){
    full_peptide<-cbind(full_peptide, t(info[match(full_peptide[,3], info[,2]),]))
   }else{
    full_peptide<-cbind(full_peptide, info[match(full_peptide[,3], info[,2]),])
   }
   colnames(full_peptide)<-tag
   
   #Get Unique Position
   uq1<-unique(info[,12])
   if(length(disable)>0){
    uq1<-uq1[-sort(match(disable,uq1))]
    full_peptide<-full_peptide[-which(!is.na(match(full_peptide[,21], disable))),]
    if(length(uq1)==0) next
   }

   #Get Unique Position with CCFP Condition
   uq2<-unique(info[as.numeric(info[,23]) > 0.8, 12])
   if(length(disable)>0){
    uq2<-uq2[-sort(match(disable,uq2))]
   }

   base_count<-nrow(min_peptide) / length(FASTA)
   unq_hla<-sapply(unique(min_peptide[,1]), 
             function(x) length(which(!is.na(match(min_peptide[,1],x))))) / base_count
   unq_hla<-unlist(sapply(1:length(unique(min_peptide[,1])), 
             function(x) rep(unique(min_peptide[,1])[x], unq_hla[x])))

   #Bind Min Peptide and Info
   if(nrow(min_peptide) > 1){
    min_peptide<-cbind(min_peptide, info[match(min_peptide[,3], info[,2]),])
    tmp<-NULL
    #HLA
    for(uqh in unq_hla){
     index<-which(!is.na(match(min_peptide[,1], uqh)))
     if(length(index)==0|length(uq1[!is.na(match(uq1, min_peptide[index,21]))])==0)next
     tmp<-rbind(tmp, min_peptide[index[sapply(uq1[!is.na(match(uq1, min_peptide[index,21]))], 
           function(x) which(!is.na(match(min_peptide[index,21], x)))
            [order(as.numeric(min_peptide[which(!is.na(match(min_peptide[index,21], x))),5]))[1]])],])
    }
    min_peptide<-tmp
   }else{
    min_peptide<-cbind(min_peptide, t(info[match(min_peptide[,3], info[,2]),]))
   }
   colnames(min_peptide)<-tag

   #Bind Min_Upper Peptide and Info
   if(is.null(min_peptide_50)) min_peptide_50<-matrix(nrow=0, ncol=length(tag) - ncol(info), 0)
   if(nrow(min_peptide_50) > 1){
    min_peptide_50<-cbind(min_peptide_50, info[match(min_peptide_50[,3], info[,2]),])
    tmp<-NULL
    for(uqh in unq_hla){
     index<-which(!is.na(match(min_peptide_50[,1], uqh)))
     if(length(index)==0|length(uq1[!is.na(match(uq1, min_peptide_50[index,21]))])==0)next
     tmp<-rbind(tmp, min_peptide_50[index[sapply(uq1[!is.na(match(uq1, min_peptide_50[index,21]))], 
           function(x) which(!is.na(match(min_peptide_50[index,21], x)))
            [order(as.numeric(min_peptide_50[which(!is.na(match(min_peptide_50[index,21], x))),5]))[1]])],])
    }
    min_peptide_50<-tmp
   } else {
    if(nrow(min_peptide_50)==0) {min_peptide_50<-cbind(min_peptide_50, info[match(min_peptide_50[,3], info[,2]),])
    } else {min_peptide_50<-cbind(min_peptide_50, t(info[match(min_peptide_50[,3], info[,2]),]))}
   }
   if(is.null(min_peptide_50)) min_peptide_50<-matrix(nrow=0, ncol=length(tag), 0)
   colnames(min_peptide_50)<-tag

   #Bind Rank Peptide and Info
   if(nrow(rank_peptide) > 1){
    rank_peptide<-cbind(rank_peptide, info[match(rank_peptide[,3], info[,2]),])
    tmp<-NULL
    for(uqh in unq_hla){
     index<-which(!is.na(match(rank_peptide[,1], uqh)))
     if(length(index)==0|length(uq1[!is.na(match(uq1, rank_peptide[index,21]))])==0)next
     tmp<-rbind(tmp, rank_peptide[index[sapply(uq1[!is.na(match(uq1, rank_peptide[index,21]))], 
           function(x) which(!is.na(match(rank_peptide[index,21], x)))
            [order(as.numeric(rank_peptide[which(!is.na(match(rank_peptide[index,21], x))),6]))[1]])],])
    }
    rank_peptide<-tmp
   } else {
    rank_peptide<-cbind(rank_peptide, t(info[match(rank_peptide[,3], info[,2]),]))
   }
   colnames(rank_peptide)<-tag

   #Bind Rank_Upper Peptide and info
   if(is.null(rank_peptide_50)) rank_peptide_50<-matrix(nrow=0, ncol=length(tag) - ncol(info), 0)
   if(nrow(rank_peptide_50) > 1){
    rank_peptide_50<-cbind(rank_peptide_50, info[match(rank_peptide_50[,3], info[,2]),])
    tmp<-NULL
    for(uqh in unq_hla){
     index<-which(!is.na(match(rank_peptide_50[,1], uqh)))
     if(length(index)==0|length(uq1[!is.na(match(uq1, rank_peptide_50[index,21]))])==0)next
     tmp<-rbind(tmp, rank_peptide_50[index[sapply(uq1[!is.na(match(uq1, rank_peptide_50[index,21]))], 
           function(x) which(!is.na(match(rank_peptide_50[index,21], x)))
            [order(as.numeric(rank_peptide_50[which(!is.na(match(rank_peptide_50[index,21], x))),6]))[1]])],])
    }
    rank_peptide_50<-tmp
   } else {
    if(nrow(rank_peptide_50)==0){rank_peptide_50<-cbind(rank_peptide_50, info[match(rank_peptide_50[,3], info[,2]),])
    } else {rank_peptide_50<-cbind(rank_peptide_50, t(info[match(rank_peptide_50[,3], info[,2]),]))}
   }
   if(is.null(rank_peptide_50)) rank_peptide_50<-matrix(nrow=0, ncol=length(tag), 0)
   colnames(rank_peptide_50)<-tag

   #Basic Information
   temp<-c(prefix,
      paste(unique(full_peptide[,1]), collapse=";"),
         scan(paste(dir2, files[grep("HLAtype.txt", files)], sep="/"), "character", sep="\n")[3],
    	    length(uq1),
    	       HLA_count)

   tag<-NULL
   ic50_col<-match(c("Mut_IC50"), colnames(full_peptide))
   rank_col<-match(c("Mut_Rank"), colnames(full_peptide))
   fpkm_col<-match(c("TumorRNA"), colnames(full_peptide))
   for(ccfp in c(0, 0.8)){
    for(key1 in c("FPKMNO", "FPKMUP", "FPKMDOWN")){
     if(key1=="FPKMNO") que<-c("")
     if(key1=="FPKMUP") que<-c(0.0001, 0.1, 1, 2, 5)
     if(key1=="FPKMDOWN") que<-c(0.0001, 0.1, 1)
     for(fpkm in que){
      #CCFP
      if(ccfp == 0) {
       ccfp_flag<-rep(TRUE, nrow(info))
      } else {
       ccfp_flag<-!is.na(as.numeric(info[,23])) & as.numeric(info[,23]) > ccfp
      }
      #FPKM
      if(key1=="FPKMNO"){
       fpkm_flag<-rep(TRUE, nrow(info))
      } else if(key1=="FPKMUP"){
       fpkm_flag<-!is.na(as.numeric(info[,19])) & as.numeric(info[,19]) > fpkm
      } else if(key1=="FPKMDOWN"){
       fpkm_flag<-!is.na(as.numeric(info[,19])) & as.numeric(info[,19]) < fpkm
      }
      
      uqx<-unique(info[which(ccfp_flag & fpkm_flag), 12])
      if(length(disable)>0){
       uqx<-uqx[-sort(match(disable,uqx))]
      }

      #Column
      tag<-c(tag, paste("Length-without", sep="-"))
      temp<-c(temp, length(uqx))

      for(ic50 in c(50, 100, 200, 500)){
       tag<-c(tag, 
           paste("IC50", ic50, key1, fpkm, "CCFP", ccfp, "Mer", sep="-"), 
	         paste("IC50", ic50, key1, fpkm, "CCFP", ccfp, "Sep", sep="-"),
           paste("IC50_wt500", ic50, key1, fpkm, "CCFP", ccfp, "Mer", sep="-"), 
	         paste("IC50_wt500", ic50, key1, fpkm, "CCFP", ccfp, "Sep", sep="-"))

       temp<-c(temp,
       #min50
       length(which(sapply(uqx, function(x) min(as.numeric(min_peptide[!is.na(match(min_peptide[,21], x)), ic50_col]))) <= ic50)),
       length(which(as.numeric(min_peptide[!is.na(match(min_peptide[,21], uqx)), ic50_col]) <= ic50)),
       #50-500
       length(which(sapply(uqx, function(x) min(as.numeric(min_peptide_50[!is.na(match(min_peptide_50[,21], x)), ic50_col]))) <= ic50)),
       length(which(as.numeric(min_peptide_50[!is.na(match(min_peptide_50[,21], uqx)), ic50_col]) <= ic50)))
      }
      
      for(rank in c(0.1, 0.5, 2)){
       tag<-c(tag, 
           paste("rank", rank, key1, fpkm, "CCFP", ccfp, "Mer", sep="-"), 
	         paste("rank", rank, key1, fpkm, "CCFP", ccfp, "Sep", sep="-"),
           paste("rank_wt2", rank, key1, fpkm, "CCFP", ccfp, "Mer", sep="-"), 
	         paste("rank_wt2", rank, key1, fpkm, "CCFP", ccfp, "Sep", sep="-"))
       temp<-c(temp,
       #min50
       length(which(sapply(uqx, function(x) min(as.numeric(rank_peptide[!is.na(match(rank_peptide[,21], x)), rank_col]))) <= rank)),
       length(which(as.numeric(rank_peptide[!is.na(match(rank_peptide[,21], uqx)), rank_col]) <= rank)),
       #50-500
       length(which(sapply(uqx, function(x) min(as.numeric(rank_peptide_50[!is.na(match(rank_peptide_50[,21], x)), rank_col]))) <= rank)),
       length(which(as.numeric(rank_peptide_50[!is.na(match(rank_peptide_50[,21], uqx)),rank_col]) <= rank)))
      }
     }
    }
   }

  e<-cbind(e, temp)
  }
  write.table(full_peptide, paste(dir2, "/", prefix, ".ALL.txt", sep=""), row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")
  write.table(min_peptide, paste(dir2, "/", prefix, ".min_lower50.txt", sep=""), row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")
  write.table(min_peptide_50, paste(dir2, "/", prefix, ".min_upper500.txt", sep=""), row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")
  write.table(rank_peptide, paste(dir2, "/", prefix, ".rank_lower0.5.txt", sep=""), row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")
  write.table(rank_peptide_50, paste(dir2, "/", prefix, ".rank_upper2.txt", sep=""), row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")

  if(is.null(e)) next
  rownames(e)<-c("HLA Type","HLA","Tumor","Nonsynonymous","Valid_HLA_Count", tag)
  write.table(e, paste(dir2, "/", prefix, ".count.txt", sep=""), 
     row.names=TRUE, col.names=FALSE, quote=FALSE, sep="\t")
}
