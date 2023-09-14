MergeFragmentsClass2<-function(hmdir = getwd(),
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
  cinfo<-c("", "Gene_ID", "Chr", "NM_ID", "ReadingFrame", "SequenceNumber", "Chrs", "NM_IDs", "GeneIDs", "Exon_Starts",
           "Exon_Ends", "GroupID", "NumOfPeptides", "NumOfStops", "Wt_Peptide", "Mutant_Peptide",
           "Wt_DNA", "Mutant_DNA", "Total_RNA", "Tumor_RNA_Ratio", "Tumor_RNA", "Tumor_RNA_based_on_DNA",
           "nB", "Checker", "MutRatio", "MutRatio_Min", "MutRatio_Max")
  info<-info[, 1:length(cinfo)]

  if(is.null(ncol(info))) info<-t(as.matrix(info))
  rownames(info)<-NULL
  colnames(info)<-cinfo

  info[,12]<-paste(info[,3], info[,12], sep="_")

  #Remove RNAseq Info
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
    test1 <- read_1col_by_fread_or_scan(paste(dir, f, sep="/"))
    test1<-gsub(" <=WB| <=SB", "", test1)
    ss1<-grep("Pos ", test1) + 2
    ee1<-grep("of strong", test1) - 2
    num1<-sapply(gsub("[ ]+", "\t", test1[ss1]), function(x) strsplit(x, "\t")[[1]][4])

    #if(length(grep("No peptides derived", test1[1:45])) > 0) next
    if(length(grep("cannot be found in hla_pseudo list", test1)) > 0) next
    if(length(grep("Could not find allele", test1)) > 0) next
    for(h1 in 1:length(num1)){
      print(paste((h1 / length(num1)) * 100, "perc. fin"))
      if(ss1[h1] == ee1[h1]){
        d1<-t(strsplit(gsub("[ ]+", "\t", test1[ss1[h1]:ee1[h1]]), "\t")[[1]][c(2, 5, 4, 6, 3, 9, 10)])
        d1<-t(d1[sapply(d1[, 5], function(x) length(grep(x, info[match(num1[h1], info[, 2]), 15]))==0),])
      } else {
        d1<-t(sapply(gsub("[ ]+", "\t", test1[ss1[h1]:ee1[h1]]), function(x) strsplit(x, "\t")[[1]][c(2, 5, 4, 6, 3, 9, 10)]))
        d1<-d1[sapply(d1[, 5], function(x) length(grep(x, info[match(num1[h1], info[, 2]), 15]))==0),]
        if(is.null(nrow(d1))) d1<-t(d1)
      }
      if(nrow(d1)==0 | ncol(d1)==0) {
        r_can<-match(num1[h1], info[,2])
        if(is.na(r_can)){r_can<-grep(num1[h1], info[,2])}
        remove<-c(remove, r_can)
        next
      }
      rownames(d1) <- NULL
      full_peptide<-rbind(full_peptide, d1)
    }
  }

  if(is.null(full_peptide)) return(NULL)
  if(nrow(full_peptide)==0) return(NULL)

  #Bind Full Peptide and info
  tag<-c("HLA", "Pos", "Gene", "Evaluated_Mutant_Peptide_Core", "Evaluated_Mutant_Peptide", "Mut_EL",
         "Mut_Rank", "Chr", "NM_ID", "ReadingFrame", "SequenceNumber", "Chrs", "NM_IDs", "GeneIDs",
         "Exon_Starts", "Exon_Ends", "GroupID", "NumOfPeptides", "NumOfStops", "Wt_Peptide",
         "Mutant_Peptide", "Total_RNA", "Tumor_RNA_Ratio", "Tumor_RNA", "Tumor_RNA_based_on_DNA",
         "MutRatio", "MutRatio_Min", "MutRatio_Max")
  colnames(full_peptide)<-tag[1:ncol(full_peptide)]
  if(nrow(full_peptide)==1){
    full_peptide<-cbind(full_peptide, t(info[match(substr(full_peptide[, 3], 1, 10), substr(info[, 2], 1, 10)),]))
  } else {
    full_peptide<-cbind(full_peptide, info[match(substr(full_peptide[, 3], 1, 10), substr(info[, 2], 1, 10)),])
  }

  full_peptide<-full_peptide[,match(tag, colnames(full_peptide))]
  write.table(full_peptide, paste(dir, "/", file_prefix, ".HLACLASS2.ALL.txt", sep=""),
              row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")
  return(full_peptide)
}
