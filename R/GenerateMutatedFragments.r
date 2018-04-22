GenerateMutatedFragments<-function(input_sequence,
                                   hmdir,
                                   job_id,
                                   refflat_file,
                                   refmrna_file,
                                   max_peptide_length,
                                   min_peptide_length,
                                   reading_frame,
                                   export_dir){

  #READ refFlat
  list_nm <- fread(refflat_file, stringsAsFactors=FALSE, sep="\n", data.table = FALSE)
  tmp <- t(apply(list_nm, 1, function(x) strsplit(x[1], "\t")[[1]]))
  list_nm_gene <- tmp[,1]
  list_nm_cut <- tmp[,2]

  #Get RNA-Code Data
  list_mra <- scan(refmrna_file, "character", sep=" ")
  list_mra <- fread(refmrna_file, stringsAsFactors=FALSE, sep='\t', data.table = FALSE)

    start_<-grep(">", list_mra)
  end_<-c(start_[-1] - 1, length(list_mra))
  list_fl_NMID<-gsub(">", "", list_mra[start_])
  list_fl_dna <-sapply(1:length(start_), function(x) paste(list_mra[(start_[x] + 2):end_[x]], collapse = ""))

  trans_from<-c("a", "t", "g", "c")
  trans_to<-c("t", "a", "c", "g")

  fasta<-NULL
  refFasta<-NULL
  random<-0

  #Obtain refFLAT Data
  s_variants <- match(nm_id, list_nm_cut)
  if(is.na(s_variants)) {
    s_variants <- which(!is.na(match(list_nm_gene, gene_symbol)))
    print(paste("NM_ID NOT Macth, Skip:", nm_id))
    next
  }

  #Calculate Sets for NM_ID, because NM_id:ExonRegion is not unique!!
  for(v in s_variants){
    #Whether Last or Not
    nm_sep <- strsplit(list_nm[v], "\t")[[1]]
    nm_id <- nm_sep[2]

    #Skip Such As "ch5_hap"
    if(nchar(nm_sep[3]) > 5) next
    strand <- nm_sep[4]
    g_name <- nm_sep[1]

    #Get Translation Start/End, Exon Start/End
    trans_start<-as.numeric(nm_sep[7])
    trans_end<-as.numeric(nm_sep[8])
    exon_start<-as.numeric(strsplit(nm_sep[10], ",")[[1]])
    exon_end<-as.numeric(strsplit(nm_sep[11], ",")[[1]])

    #Obtain DNA sequence of Transcriptome
    #DNAseq is Unique
    dna<-list_fl_dna[match(nm_id, list_fl_NMID)]
    if(nchar(dna) != sum(exon_end - exon_start)){
      dif <- sum(exon_end - exon_start) - nchar(dna)
      if(abs(dif) <= 0) {
        print(paste("cDNA Length does not Match to Exon-Start/End Length, Skip", nm_id))
        next
      }
    }

    #Get Relative Translation-Start Position (0-start to 1-start)
    if(strand=="+"){
      point<-(exon_end > trans_start)
      ts_point<-sum((exon_end - exon_start)[!point]) +
        (trans_start - exon_start[which(point)[1]]) + 1
    }else{
      point<-(exon_start > trans_end)
      ts_point<-sum((exon_end - exon_start)[point]) +
        (exon_end[rev(which(!point))[1]] - trans_end) + 1
    }

    #Check Start Codon
    if(substr(dna, ts_point, ts_point + 2)!="atg"){
      print(paste("Start Position is not ATG, Skip", nm_id))
      next
    }

    #Get Relative Translation-End Position
    if(strand=="+"){
      point<-(exon_end >= trans_end)
      te_point<-sum((exon_end - exon_start)[!point]) +
        (trans_end - exon_start[which(point)[1]])
    } else {
      point<-(exon_start > trans_start)
      te_point<-sum((exon_end - exon_start)[point]) +
        (exon_end[rev(which(!point))[1]] - trans_start)
    }

    #Check Stop Codon
    if(amino[match(substr(dna, te_point-2, te_point), codon)]!="X"){
        print(paste("End Position Amino Acid is not X, Skip", nm_id))
        next
    }

    #Check Peptide Length
    stop_loop<-FALSE
    dna_trans <- substr(dna, ts_point, te_point)

    #Translation Region is not Valid
    if(nchar(dna_trans)%%3!=0) {
      print("The Length of RNA is not a multiple of 3, Skip")
      next
    }

    #Make Normal Peptide
    peptide_normal <- NULL
    dna_trans_normal <- dna_trans
    while(nchar(dna_trans)>=3){
      peptide_normal <- c(peptide_normal,amino[match(substr(dna_trans, 1, 3), codon)])
      dna_trans <- substr(dna_trans, 4, nchar(dna_trans))
    }
    if(match("X", peptide_normal) < length(peptide_normal)){
      next
    }

    #Make Mutated Peptide
    peptide_mutated <- NULL
    input_sequence <- substr(tolower(input_sequence), reading_frame, nchar(input_sequence))
    input_sequence_2 <- input_sequence
    while(nchar(input_sequence_2)>=3){
      peptide_mutated <- c(peptide_mutated, amino[match(substr(input_sequence_2, 1, 3), codon)])
      input_sequence_2 <- substr(input_sequence_2, 4, nchar(input_sequence_2))
    }
    if(!is.na(match("X", peptide_mutated))){
      if(match("X", peptide_mutated) < length(peptide_mutated)){
        peptide_mutated <- peptide_mutated[1:(match("X", peptide_mutated) - 1)]
      }
    }
    pep_end_pos <- length(peptide_mutated) - min_peptide_length
    if(pep_end_pos < 1) {
      print("Mutated Peptide is too short, Skip")
      next
    }
    flg_vec <- rep(FALSE, length(peptide_mutated))
    ref_pep <- paste(peptide_normal, collapse = "")
    for(i in 1:(pep_end_pos + 1)){
      flg <- length(grep(paste(peptide_mutated[i:(i + min_peptide_length - 1)], collapse = ""), ref_pep)) == 0
      if(flg) flg_vec[(ifelse(i - (max_peptide_length - min_peptide_length) < 1, 1, i - (max_peptide_length - min_peptide_length))):
                       (ifelse(i + max_peptide_length - 1 > length(peptide_mutated), length(peptide_mutated), i + max_peptide_length - 1))] <- TRUE
    }

    peptide_mutated <- paste(ifelse(flg_vec, peptide_mutated, "-"), collapse = "")
    for(peptide in strsplit(peptide_mutated, "-")[[1]]){
      #Save Peptide
      if(nchar(peptide) < min_peptide_length) break
      refFasta<-rbind(refFasta,
                      c(paste(random, gsub("\"","", g_name), sep="_"),
                        nm_sep[3],
                        nm_id,
                        NA,
                        NA,
                        NA,
                        NA,
                        NA,
                        exon_start[1],
                        rev(exon_end)[1],
                        NA,
                        NA,
                        NA,
                        paste(peptide_normal, collapse=""),
                        peptide,
                        dna_trans_normal,
                        input_sequence))

        #Remove X and Save Fasta in Mutated Peptide
        if(!is.na(match("X", peptide))){
          peptide <- peptide[1:(match("X", peptide) - 1)]
        }
        fasta <- c(fasta, sub("_","", paste(">", random, gsub("\"","", g_name), sep="_")))
        fasta <- c(fasta, paste(peptide, collapse=""))

        random <- random + 1
        print("Peptide Successfully Generated!!")
    }
  }

  #Integrate The Same Peptide
  if(is.null(refFasta)) {
    return(NULL)
  }
  i<-1
  if(!is.null(nrow(refFasta))){
    while(i<=nrow(refFasta)){
      #If mutation position, mutated peptide, normal peptide are all the same, Integrate These Peptides
      #Note That, If the mutation position and chromosome number are the same, Merge script integrate them in the later process.
      hit <- which(refFasta[i, 14]==refFasta[,14] &
                   refFasta[i, 15]==refFasta[,15])
      if(length(hit)==1){
        i<-i+1
        next
      }
      #Collapse NM_ID, Info, Exon-start, Exon-end
      temp1<-paste(refFasta[hit,3],collapse=";")
      temp2<-paste(refFasta[hit,4],collapse=";")
      temp3<-paste(refFasta[hit,9],collapse=";")
      temp4<-paste(refFasta[hit,10],collapse=";")
      refFasta[i,3]<-temp1
      refFasta[i,4]<-temp2
      refFasta[i,9]<-temp3
      refFasta[i,10]<-temp4
      refFasta<-refFasta[-hit[-1],]
      if(is.null(nrow(refFasta))){
        refFasta<-t(refFasta)
      }
      fasta<-fasta[-(c(hit[-1] * 2 - 1, hit[-1] * 2))]
      i<-i+1
    }
  }
  write.table(fasta,
              paste(export_dir, "/", job_id, ".", "peptide", ".", "fasta", sep=""),
              row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")
  write.table(cbind(refFasta, matrix(nrow = nrow(refFasta), ncol = 27 - ncol(refFasta), NA)),
              paste(export_dir, "/", job_id, ".", "peptide", ".", "txt", sep=""),
              row.names=seq(1:nrow(refFasta)), col.names=FALSE, quote=FALSE, sep="\t")
}
