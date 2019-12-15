GenerateMutatedFragments<-function(input_sequence,
                                   input_nm_id,
                                   group_ids,
                                   hmdir,
                                   job_id,
                                   refflat_file,
                                   refmrna_file,
                                   max_peptide_length,
                                   min_peptide_length,
                                   reading_frame,
                                   export_dir,
                                   reference_nm_id,
                                   reference_gene_symbol,
                                   IgnoreShortPeptides){

  #READ refFlat
  list_nm <- read_refFlat(refflat_file)
  list_nm_gene <- list_nm[, 1]
  list_nm_cut <- list_nm[, 2]

  #Get RNA-Code Data
  list_mra <- read_refmrn(refmrna_file)
  start_ <- grep(">", list_mra)
  end_ <- c(start_[-1] - 1, length(list_mra))
  list_fl_NMID <- gsub(">", "", sapply(list_mra[start_], function(x) strsplit(x, " ")[[1]][1]))
  list_fl_dna <- sapply(1:length(start_), function(x) paste(list_mra[(start_[x] + 1):end_[x]], collapse = ""))

  trans_from<-c("a", "t", "g", "c")
  trans_to<-c("t", "a", "c", "g")

  #Obtain refFLAT Data
  s_variants_from_nmid <- match(reference_nm_id, list_nm_cut)
  s_variants_from_gene <- which(!is.na(match(list_nm_gene, reference_gene_symbol)))
  s_variants <- sort(unique(c(s_variants_from_nmid, s_variants_from_gene)))
  if(length(s_variants) == 0){
    print("No Wt-NM_ID Identified")
  }

  #Calculate Sets for NM_ID, because NM_id:ExonRegion is not unique!!
  peptide_normal_merged <- NULL
  dna_trans_normal_merged <- NULL
  chrs <- NULL
  gene_ids <- NULL
  nm_ids <- NULL
  exon_starts <- NULL
  exon_ends <- NULL
  for(v in s_variants){
    #Whether Last or Not
    nm_sep <- sapply(list_nm[v, ], as.character)
    nm_id <- nm_sep[2]
    nm_ids <- paste(nm_ids, nm_id, sep = ifelse(length(nm_ids) > 0, ";", ""))

    #Skip Such As "ch5_hap"
    if(nchar(nm_sep[3]) > 5) next
    chr <- nm_sep[3]
    chrs <- paste(chrs, nm_sep[3], sep = ifelse(length(chrs) > 0, ";", ""))
    strand <- nm_sep[4]
    g_name <- nm_sep[1]
    gene_ids <- paste(gene_ids, g_name, sep = ifelse(length(gene_ids) > 0, ";", ""))

    #Get Translation Start/End, Exon Start/End
    trans_start <- as.numeric(nm_sep[7])
    trans_end <- as.numeric(nm_sep[8])
    exon_start <- as.numeric(strsplit(nm_sep[10], ",")[[1]])
    exon_starts <- paste(exon_starts, exon_start[1], sep = ifelse(length(exon_starts) > 0, ";", ""))
    exon_end <- as.numeric(strsplit(nm_sep[11], ",")[[1]])
    exon_ends <- paste(exon_ends, rev(exon_end)[1], sep = ifelse(length(exon_ends) > 0, ";", ""))

    #Obtain DNA sequence of Transcriptome
    dna <- list_fl_dna[match(nm_id, list_fl_NMID)]

    #Check DNA
    if(check_dna_validity(dna, nm_id, exon_end, exon_start, ambiguous_between_exon, final_s_variants, Pass)) next

    #Get Relative Translation-Start Position (0-start to 1-start)
    ts_point <- get_relative_translation_start_position(strand, exon_end, trans_start, exon_start, trans_end)

    #Check Start Codon
    d <- check_start_codon(dna, ts_point, ambiguous_codon, nm_id)
    if(d < -998 | is.null(d)) next

    #Get Relative Translation-End Position
    te_point <- get_relative_translation_end_position(strand, exon_end, trans_start, exon_start, trans_end)

    #Check Stop Codon
    e <- check_stop_codon(dna, te_point, ts_point, ambiguous_codon, amino, nm_id)
    if(e < -998 | is.null(e)) next

    #Check Peptide Length
    stop_loop<-FALSE
    dna_trans <- substr(dna, ts_point, te_point)

    #Translation Region is not Valid
    if(nchar(dna_trans)%%3!=0) {
      print("Translation Region is not Valid.")
      next
    }

    #Make Normal Peptide
    dna_trans_normal_merged <- c(dna_trans_normal_merged, dna_trans)
    peptide_normal <- make_normal_peptide(dna_trans, amino, codon, 0, 0)
    if(is.null(peptide_normal)) next
    peptide_normal_merged <- c(peptide_normal_merged, ifelse(length(peptide_normal_merged) > 0, "X", ""), peptide_normal)
  }

  #Make Mutated Peptide
  peptide_mutated <- NULL

  #Obtain refFLAT Data
  s_variants_from_input_nmid <- NULL
  if(!is.na(input_nm_id[1])){
    s_variants_from_input_nmid <- match(input_nm_id, list_nm_cut)
    s_variants_from_input_nmid <- sort(unique(s_variants_from_input_nmid))
  }
  if(is.na(input_sequence[1]) & length(s_variants_from_input_nmid) == 0){
    print(paste("NM_ID NOT Macth, Skip:", nm_id))
    return(NULL)
  }

  #
  names(input_sequence) <- rep("", length(input_sequence))
  if(length(s_variants_from_input_nmid) != 0){
    if(is.na(input_sequence[1])) input_sequence <- NULL
    for(v in s_variants_from_input_nmid){
      #Whether Last or Not
      nm_sep <- as.character(list_nm[v,])
      nm_id <- nm_sep[2]

      #Skip Such As "ch5_hap"
      if(nchar(nm_sep[3]) > 5) next
      chr <- nm_sep[3]
      strand <- nm_sep[4]
      g_name <- nm_sep[1]

      #Get Translation Start/End, Exon Start/End
      trans_start <- as.numeric(nm_sep[7])
      trans_end <- as.numeric(nm_sep[8])
      exon_start <- as.numeric(strsplit(nm_sep[10], ",")[[1]])
      exon_end <- as.numeric(strsplit(nm_sep[11], ",")[[1]])

      #Obtain DNA sequence of Transcriptome
      dna <- list_fl_dna[match(nm_id, list_fl_NMID)]

      #Check DNA
      if(check_dna_validity(dna, nm_id, exon_end, exon_start, ambiguous_between_exon, final_s_variants, Pass)) next

      #Get Relative Translation-Start Position (0-start to 1-start)
      ts_point <- get_relative_translation_start_position(strand, exon_end, trans_start, exon_start, trans_end)

      #Check Start Codon
      d <- check_start_codon(dna, ts_point, ambiguous_codon, nm_id)
      if(d < -998 | is.null(d)) next

      #Get Relative Translation-End Position
      te_point <- get_relative_translation_end_position(strand, exon_end, trans_start, exon_start, trans_end)

      #Check Stop Codon
      e <- check_stop_codon(dna, te_point, ts_point, ambiguous_codon, amino, nm_id)
      if(e < -998 | is.null(e)) next

      #Check Peptide Length
      stop_loop<-FALSE
      dna_trans <- substr(dna, ts_point, nchar(dna))

      count_dna <- 0
      dna_trans <- substr(tolower(dna_trans), reading_frame, nchar(dna_trans))
      while(nchar(dna_trans) >= 3){
        if(amino[match(substr(dna_trans, 1, 3), codon)] == "X" & count_dna >= te_point - ts_point - 3) break
        dna_trans <- substr(dna_trans, 4, nchar(dna_trans))
        count_dna <- count_dna + 3
      }
      dna_trans <- substr(dna, ts_point, min(ts_point + count_dna + 3, nchar(dna)))

      input_sequence <- c(input_sequence, dna_trans)
      names(input_sequence)[length(input_sequence)] <- paste(g_name, nm_id, chr, sep = ";")
    }
  }

  if(is.null(input_sequence[1])) return(NULL)
  if(is.na(input_sequence[1])) return(NULL)

  fasta<-NULL
  refFasta<-NULL
  random<-0
  for(input_sequence_1 in input_sequence){
    peptide_mutated <- NULL
    input_sequence_2 <- substr(tolower(input_sequence_1), reading_frame, nchar(input_sequence_1))
    while(nchar(input_sequence_2) >= 3){
      peptide_mutated <- c(peptide_mutated, amino[match(substr(input_sequence_2, 1, 3), codon)])
      input_sequence_2 <- substr(input_sequence_2, 4, nchar(input_sequence_2))
    }

    number_of_peptide <- length(peptide_mutated)
    peptide_mutated_sep <- strsplit(paste(peptide_mutated, collapse = ""), "X")[[1]]
    number_of_stop <- length(peptide_mutated_sep) - 1
    for(peptide_mutated in peptide_mutated_sep){
      peptide <- peptide_mutated
      peptide_mutated <- strsplit(peptide_mutated, "")[[1]]
      # if(length(peptide_mutated) >= min_peptide_length) {
      #   pep_end_pos <- length(peptide_mutated) - min_peptide_length
      #   if(pep_end_pos < 1) {
      #     print("Mutated Peptide is too short, Skip")
      #     next
      #   }
      #   flg_vec <- rep(FALSE, length(peptide_mutated))
      #   ref_pep <- paste(peptide_normal_merged, collapse = "")
      #   for(i in 1:(pep_end_pos + 1)){
      #     flg <- length(grep(paste(peptide_mutated[i:(i + min_peptide_length - 1)], collapse = ""), ref_pep)) == 0
      #     if(flg) flg_vec[(ifelse(i - (max_peptide_length - min_peptide_length) < 1, 1, i - (max_peptide_length - min_peptide_length))):
      #                       (ifelse(i + max_peptide_length - 1 > length(peptide_mutated), length(peptide_mutated), i + max_peptide_length - 1))] <- TRUE
      #   }
      #   peptide_mutated <- paste(ifelse(flg_vec, peptide_mutated, "-"), collapse = "")
      # } else {
      #   peptide_mutated <- paste(peptide_mutated, collapse = "")
      # }
      # for(peptide in strsplit(peptide_mutated, "-")[[1]]){
        #Save Peptide
        #if(nchar(peptide) < min_peptide_length) next
        seq_num <- match(input_sequence_1, input_sequence)
        g_name <- strsplit(names(input_sequence)[seq_num], ";")[[1]][1]
        g_name <- ifelse(is.na(g_name) | g_name == "", substr(input_sequence_1, 1, 10), g_name)
        nm_id <- strsplit(names(input_sequence)[seq_num], ";")[[1]][2]
        nm_id <- ifelse(is.na(nm_id), "", nm_id)
        group_id <- group_ids[match(input_sequence_1, names(group_ids))]
        if(is.na(group_id)) group_id <- group_ids[match(nm_id, names(group_ids))]

        refFasta<-rbind(refFasta,
                        c(paste(random, g_name, sep="_"),
                        0,
                        nm_id,
                        reading_frame,
                        seq_num,
                        ifelse(is.null(chrs), "", chrs),
                        ifelse(is.null(nm_ids), "", nm_ids),
                        ifelse(is.null(gene_ids), "", gene_ids),
                        ifelse(is.null(exon_starts), "", exon_starts),
                        ifelse(is.null(exon_ends), "", exon_ends),
                        group_id,
                        number_of_peptide,
                        number_of_stop,
                        ifelse(is.null(paste(peptide_normal_merged, collapse="")), "", paste(peptide_normal_merged, collapse="")),
                        paste(peptide_mutated, collapse = ""),
                        ifelse(is.null(dna_trans_normal_merged), "", dna_trans_normal_merged),
                        input_sequence_1))

        fasta <- c(fasta, sub("_","", paste(">", random, gsub("\"","", g_name), sep="_")))
        fasta <- c(fasta, paste(peptide, collapse=""))

        random <- random + 1
        print("Peptide Successfully Generated!!")

    }
  }

  write.table(fasta,
              paste(export_dir, "/", job_id, ".", "peptide", ".", "fasta", sep=""),
              row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")
  write.table(cbind(refFasta, matrix(nrow = nrow(refFasta), ncol = 27 - ncol(refFasta) - 1, NA)),
              paste(export_dir, "/", job_id, ".", "peptide", ".", "txt", sep=""),
              row.names=seq(1:nrow(refFasta)), col.names=FALSE, quote=FALSE, sep="\t")
}
