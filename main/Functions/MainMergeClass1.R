MainMergeClass1<-function(hmdir = getwd(), input_dir, input_file_prefix, Tumor_RNA_BASED_ON_DNA = TRUE, INDEL = FALSE){
  dir<-paste(hmdir, input_dir, sep="/")
  files<-list.files(paste(dir, sep="/"))
  if(INDEL){
    files<-files[grep("INDEL", files)]
  } else {
    files<-files[grep("INDEL", files, invert = TRUE)]
  }
  
  #Get Peptide Info
  files_part<-files[intersect(grep("HLACLASS1", files), grep(input_file_prefix, files))]
  if(length(files_part)==0){
    print("No File:!!")
    return(NULL)
  }
  files_rest<-files[grep("HLACLASS", files, invert = TRUE)]
  file_peptide_info<-files_rest[intersect(grep(input_file_prefix, files_rest), grep("peptide\\.txt$", files_rest))]
  info<-t(sapply(scan(paste(dir, file_peptide_info, sep="/"), "character", sep="\n"), function(x) strsplit(x, "\t")[[1]]))
  cinfo<-c("", "Gene ID", "Chr", "NM_ID", "Change", "ref", "alt", "Prob", "Mutation Prob.", "Exon Start", 
           "Exon End", "Mutation Position", "Depth", "TumorDepth", "Peptide Normal", "Peptide Mutation", 
           "DNA_Normal", "DNA_Mut", "TotalRNA", "TumorRNARatio", "TumorRNA", "nA", "nB", "Checker", 
           "MutRatio", "MutRatio Min", "MutRatio Max")
  info<-info[,1:length(cinfo)]
  
  if(is.null(ncol(info))){info<-t(as.matrix(info))}
  row.names(info)<-NULL
  colnames(info)<-cinfo
  
  if(Tumor_RNA_BASED_ON_DNA){
    info[,match("TumorRNA",colnames(info))]<-as.numeric(info[,match("TotalRNA",colnames(info))]) * as.numeric(info[,match("TumorDepth",colnames(info))]) /
      as.numeric(info[,match("Depth",colnames(info))])
  }
  
  #Remove RNAseq Info
  rownames(info)<-NULL
  info<-info[,-match(c("DNA_Normal", "DNA_Mut"), colnames(info))]
  if(is.null(ncol(info))){info<-t(as.matrix(info))}
  
  #Include Stop Codon
  remove<-which(sapply(info[,c(16)], function(x) length(grep("X", rev(strsplit(x, "")[[1]])[-1]))>0))
  if(length(remove) > 0) info<-info[-remove,]
  if(is.null(ncol(info))){info<-t(as.matrix(info))}
  if(nrow(info)==0) return(NULL)
  
  #Each,All,DRB
  full_peptide<-NULL
  min_peptide<-NULL
  min_peptide_50<-NULL
  rank_peptide<-NULL
  rank_peptide_50<-NULL
  
  #allele,start,end,length,peptide,ic50,Rank,Peptide_Normal_Sep,norm_ic_50,norm_Rank
  for(f in files_part[grep("normpeptide", files_part, invert = TRUE)]){
    print(paste(dir, f, sep="/"))
    
    test1<-scan(paste(dir, f, sep="/"), "character", sep="\n", skip=1)
    test1<-gsub(" <= WB", "", test1)
    ss1<-grep(" Pos ", test1)+2
    ee1<-grep("Protein", test1)-2
    num1<-sapply(gsub("[ ]+","\t",test1[ss1]), function(x) strsplit(x, "\t")[[1]][12])

    test2<-scan(paste(dir, sub("peptide\\.txt", "normpeptide\\.txt", f), sep="/"),"character", sep="\n",skip=1)
    test2<-gsub(" <= WB", "", test2)
    ss2<-grep(" Pos ", test2)+2
    ee2<-grep("Protein", test2)-2
    num2<-sapply(gsub("[ ]+","\t",test2[ss2]), function(x) strsplit(x, "\t")[[1]][12])

    if(length(grep("cannot be found in hla_pseudo list", test1))>0) next
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
    
  if(is.null(full_peptide)) return(NULL)
  if(nrow(full_peptide)==0) return(NULL)
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
  uq1<-unique(apply(info[,c(3,12)], 1, paste, collapse="_"))
  
  #Do not use nrow(min_peptide_50) because peptide of which IC50 > 500 are removed
  base_count<-nrow(min_peptide) / length(files_part[grep("normpeptide", files_part, invert = TRUE)])
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
      label_index<-apply2(min_peptide[index,c(12,21)], 1, function(x){paste(x, collapse="_")})
      if(length(index)==0 | length(uq1[!is.na(match(uq1, label_index))])==0)next
      tmp<-rbind(tmp, min_peptide[index[sapply(uq1[!is.na(match(uq1, label_index))], 
                                               function(x) which(!is.na(match(label_index, x)))
                                               [order(as.numeric(min_peptide[which(!is.na(match(label_index, x))),5]))[1]])],])
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
      label_index<-apply2(min_peptide_50[index,c(12,21)], 1, function(x){paste(x, collapse="_")})
      if(length(index)==0|length(uq1[!is.na(match(uq1, label_index))])==0)next
      tmp<-rbind(tmp, min_peptide_50[index[sapply(uq1[!is.na(match(uq1, label_index))], 
                                                  function(x) which(!is.na(match(label_index, x)))
                                                  [order(as.numeric(min_peptide_50[which(!is.na(match(label_index, x))),5]))[1]])],])
    }
    min_peptide_50<-tmp
  }else{
    if(nrow(min_peptide_50)==0){min_peptide_50<-cbind(min_peptide_50, info[match(min_peptide_50[,3], info[,2]),])
    }else{min_peptide_50<-cbind(min_peptide_50, t(info[match(min_peptide_50[,3], info[,2]),]))}
  }
  if(is.null(min_peptide_50)) min_peptide_50<-matrix(nrow=0, ncol=length(tag), 0)
  colnames(min_peptide_50)<-tag
  
  #Bind Rank Peptide and Info
  if(nrow(rank_peptide) > 1){
    rank_peptide<-cbind(rank_peptide, info[match(rank_peptide[,3], info[,2]),])
    tmp<-NULL
    for(uqh in unq_hla){
      index<-which(!is.na(match(rank_peptide[,1], uqh)))
      label_index<-apply2(rank_peptide[index,c(12,21)], 1, function(x){paste(x, collapse="_")})
      if(length(index)==0|length(uq1[!is.na(match(uq1, label_index))])==0)next
      tmp<-rbind(tmp, rank_peptide[index[sapply(uq1[!is.na(match(uq1, label_index))], 
                                                function(x) which(!is.na(match(label_index, x)))
                                                [order(as.numeric(rank_peptide[which(!is.na(match(label_index, x))),6]))[1]])],])
    }
    rank_peptide<-tmp
  } else{
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
      label_index<-apply2(rank_peptide_50[index,c(12,21)], 1, function(x){paste(x, collapse="_")})
      if(length(index)==0|length(uq1[!is.na(match(uq1, label_index))])==0)next
      tmp<-rbind(tmp, rank_peptide_50[index[sapply(uq1[!is.na(match(uq1, label_index))], 
                                                   function(x) which(!is.na(match(label_index, x)))
                                                   [order(as.numeric(rank_peptide_50[which(!is.na(match(label_index, x))),6]))[1]])],])
    }
    rank_peptide_50<-tmp
  } else {
    if(nrow(rank_peptide_50)==0){rank_peptide_50<-cbind(rank_peptide_50, info[match(rank_peptide_50[,3], info[,2]),])
    }else{rank_peptide_50<-cbind(rank_peptide_50, t(info[match(rank_peptide_50[,3], info[,2]),]))}
  }
  if(is.null(rank_peptide_50)) rank_peptide_50<-matrix(nrow=0, ncol=length(tag), 0)
  colnames(rank_peptide_50)<-tag
  
  write.table(full_peptide, paste(dir, "/", input_file_prefix, ".CLASS1.ALL.txt", sep=""), row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")
  write.table(min_peptide, paste(dir, "/", input_file_prefix, ".CLASS1.IC50min.txt", sep=""), row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")
  write.table(min_peptide_50, paste(dir, "/", input_file_prefix, ".CLASS1.IC50min_WTupper500.txt", sep=""), row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")
  write.table(rank_peptide, paste(dir, "/", input_file_prefix, ".CLASS1.RANKmin.txt", sep=""), row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")
  write.table(rank_peptide_50, paste(dir, "/", input_file_prefix, ".CLASS1.RANKmin_WTupper2.txt", sep=""), row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")
}
