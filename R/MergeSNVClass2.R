MergeSNVClass2<-function(hmdir = getwd(),
                         annotation_file,
                         input_dir,
                         file_prefix){
  print("Merging Results...")

  dir<-paste(hmdir, input_dir, sep="/")
  files<-list.files(paste(dir, sep="/"))

  #Get Peptide Info
  files_part<-files[intersect(grep("HLACLASS2", files), grep(file_prefix, files))]
  if(length(files_part)==0){
    print("No File Detected!!")
    return(NULL)
  }

  info<-t(sapply(scan(paste(annotation_file, sep="/"), "character", sep="\n"), function(x) strsplit(x, "\t")[[1]]))
  cinfo<-c("", "Gene_ID", "Chr", "NM_ID", "Change", "Ref", "Alt", "Prob", "Mutation_Prob.", "Exon_Start",
           "Exon_End", "Mutation_Position", "Total_Depth", "Tumor_Depth", "Wt_Peptide", "Mutant_Peptide",
           "Wt_DNA", "Mutant_DNA", "Total_RNA", "Tumor_RNA_Ratio", "Tumor_RNA", "Tumor_RNA_based_on_DNA",
           "nB", "Checker", "MutRatio", "MutRatio_Min", "MutRatio_Max")
  info<-info[, 1:length(cinfo)]

  if(is.null(ncol(info))) info<-t(as.matrix(info))
  row.names(info)<-NULL
  colnames(info)<-cinfo

  info[,12]<-paste(info[,3], info[,12], sep="_")
  info[, match("Tumor_RNA_based_on_DNA",colnames(info))]<-
    as.numeric(info[,match("Total_RNA",colnames(info))]) *
     as.numeric(info[,match("Tumor_Depth",colnames(info))]) /
      as.numeric(info[,match("Total_Depth",colnames(info))])

  #Remove RNAseq Info
  rownames(info)<-NULL
  info<-info[, -match(c("Wt_DNA", "Mutant_DNA"), colnames(info))]
  if(is.null(ncol(info))){info<-t(as.matrix(info))}

  #Include Stop Codon
  removeX<-which(sapply(info[,c(16)], function(x) length(grep("X", rev(strsplit(x, "")[[1]])[-1]))>0))
  if(length(removeX) > 0) info<-info[-remove,]
  if(is.null(ncol(info))){info<-t(as.matrix(info))}
  if(nrow(info)==0) return(NULL)

  #allele,start,end,length,peptide,ic50,Rank,Peptide_Normal_Sep,norm_ic_50,norm_Rank
  full_peptide<-NULL
  for(f in files_part[grep("\\.peptide\\.txt", files_part)]){
    print(paste(dir, f, sep="/"))

    test1<-scan(paste(dir, f, sep="/"), "character", sep="\n", skip=1)
    test1<-gsub(" <=WB| <=SB", "", test1)
    ss1<-grep(" Pos ", test1) + 2
    ee1<-grep("of strong", test1) - 2
    num1<-sapply(gsub("[ ]+","\t",test1[ss1]), function(x) strsplit(x, "\t")[[1]][5])

    test2<-scan(paste(dir, sub("peptide\\.txt", "wtpeptide\\.txt", f), sep="/"),"character", sep="\n",skip=1)
    test2<-gsub(" <=WB| <=SB", "", test2)
    ss2<-grep(" Pos ", test2) + 2
    ee2<-grep("of strong", test2) - 2
    num2<-sapply(gsub("[ ]+","\t",test2[ss2]), function(x) strsplit(x, "\t")[[1]][5])

    #if(length(grep("No peptides derived", test1[1:45]))>0) next
    if(length(grep("cannot be found in hla_pseudo list", test1))>0) next
    if(length(grep("Could not find allele", test1))>0) next
    for(h1 in 1:length(num1)){
      print(paste((h1 / length(num1)) * 100, "perc. fin"))
      #Skip if not match
      if(is.na(grep(num1[h1], info[,2])[1]))next

      d4<-NULL
      hit<-match(num1[h1], num2)
      d1<-t(sapply(gsub("[ ]+", "\t", test1[ss1[h1]:ee1[h1]]), function(x) strsplit(x, "\t")[[1]][c(2,3,4,5,10,11)]))
      d2<-t(sapply(gsub("[ ]+", "\t", test2[ss2[hit]:ee2[hit]]), function(x) strsplit(x, "\t")[[1]][c(2,3,4,5,10,11)]))
      l1<-sapply(d1[,3], nchar)
      l2<-sapply(d2[,3], nchar)
      for(r1 in unique(l1)){
        hit1<-which(l1 == r1)
        hit2<-which(l2 == r1)
        if(length(hit1) == 0) next
        if(length(hit1) == 1){
          d3<-t(c(d1[hit1,c(2,1,4,3,5,6)], d2[hit2[match(d1[hit1,1], d2[hit2,1])], c(3,5,6)]))
        } else {
          d3<-cbind(d1[hit1,c(2,1,4,3,5,6)], d2[hit2[match(d1[hit1,1], d2[hit2,1])], c(3,5,6)])
        }
        d3<-d3[d3[, 4] != d3[, 7],]
        d4<-rbind(d4, d3)
      }
      if(nrow(d4)==0) {
        print("Warning!! d4 is zero!!")
        next
      }
      full_peptide<-rbind(full_peptide, d4)
    }
  }

  if(is.null(full_peptide)) return(NULL)
  if(nrow(full_peptide)==0) return(NULL)

  #Bind Full Peptide and info
  tag<-c("HLA", "Pos", "Gene", "Evaluated_Mutant_Peptide", "Mut_IC50", "Mut_Rank",
         "Evaluated_Wt_Peptide", "Wt_IC50", "Wt_Rank", "Chr", "NM_ID", "Change",
         "Ref", "Alt", "Prob", "Mutation_Prob.", "Exon_Start", "Exon_End",
         "Mutation_Position", "Total_Depth", "Tumor_Depth", "Wt_Peptide",
         "Mutant_Peptide", "Total_RNA", "Tumor_RNA_Ratio", "Tumor_RNA",
         "Tumor_RNA_based_on_DNA", "MutRatio", "MutRatio_Min", "MutRatio_Max")
  colnames(full_peptide)<-tag[1:ncol(full_peptide)]
  if(nrow(full_peptide)==1){
    full_peptide<-cbind(full_peptide, t(info[match(substr(full_peptide[,3], 1, 10), substr(info[,2], 1, 10)),]))
  }else{
    full_peptide<-cbind(full_peptide, info[match(substr(full_peptide[,3], 1, 10), substr(info[,2], 1, 10)),])
  }

  full_peptide<-full_peptide[,match(tag, colnames(full_peptide))]
  write.table(full_peptide, paste(dir, "/", file_prefix, ".HLACLASS2.ALL.txt", sep=""),
              row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")
  return(full_peptide)
}
