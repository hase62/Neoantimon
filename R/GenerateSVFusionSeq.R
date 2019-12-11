GenerateSVFusionSeq<-function(input_file,
                             hmdir,
                             job_id,
                             refflat_file,
                             refmrna_file,
                             max_peptide_length,
                             chr_column,
                             mutation_start_column,
                             mutation_end_column,
                             mutation_ref_column,
                             mutation_alt_bnd_column,
                             nm_id_column,
                             gene_symbol_column,
                             depth_normal_column,
                             depth_tumor_column,
                             mate_id_column,
                             ambiguous_between_exon,
                             ambiguous_codon,
                             export_dir,
                             IgnoreShortPeptides){

  index<-strsplit(scan(input_file, "character", sep="\n", nlines=1), "\t")[[1]]
  if(requireNamespace("data.table", quietly=TRUE)) {
    data <- fread(input_file, stringsAsFactors=FALSE, header = TRUE, sep="\n", data.table = FALSE)[, 1]
  } else {
    data  <- scan(input_file, "character", sep = "\n", skip = 1)
  }

  if(length(data)<1){
    q("no")
  }

  mateIDs<-sapply(sapply(data, function(x) strsplit(x, "\t")[[1]][mate_id_column]),
                  function(x) strsplit(x, "_")[[1]][1])
  uIDs<-unique(mateIDs)[sapply(unique(mateIDs), function(x) length(which(!is.na(match(mateIDs, x)))) > 1)]

  #READ refFlat
  list_nm <- read_refFlat(refflat_file)
  list_nm_gene <- tmp[, 1]
  list_nm_cut <- tmp[, 2]

  #Get RNA-Code Data
  list_mra <- read_refmrn(refmrna_file)
  start_ <- grep(">", list_mra)
  end_ <- c(start_[-1] - 1, length(list_mra))
  tmp <- gsub(">", "", sapply(list_mra[start_], function(x) strsplit(x, " ")[[1]][1]))
  list_fl_NMID <- tmp
  list_fl_dna <- sapply(1:length(start_),
                        function(x) paste(list_mra[(start_[x] + 1):end_[x]], collapse = ""))

  trans_from<-c("a", "t", "g", "c")
  trans_to<-c("t", "a", "c", "g")

  random<-0
  fasta<-NULL
  refFasta<-NULL
  for(uID in uIDs){
    print(uID)
    sets<-NULL
    for(i in which(!is.na(match(mateIDs, uID)))){
      set<-NULL
      gc();gc();
      #Extract i-th Data
      #f is a row of mutation data
      #Extract i-th Data
      f<-strsplit(data[i], "\t")[[1]]

      #Chromosome
      chr<-f[chr_column]

      whether_exon_intron<-ifelse(length(grep("exon", f)) > 0,
                                  ifelse(length(grep("intron", f)) > 0,
                                         ifelse(grep("exon", f)[1] > grep("intron", f)[1], "exon", "intron"),
                                         "exon"),
                                  "intron")

      #MP:Somatic Mutation Probability
      MP<-0
      if(length(grep("MP=", f))>0){
        MP<-as.numeric(strsplit(strsplit(f[grep("MP=",f)], "MP=")[[1]][2],";")[[1]][1])
      }

      #GP:Genotype Probability
      GP<-0
      if(length(grep("GP=",f))>0){
        GP<-strsplit(strsplit(f[grep("GP=",f)], "GP=")[[1]][2],";")[[1]][1]
      }

      #DP:Total Depth
      DP<-0
      alt<-NULL
      ref<-NULL
      if(!is.na(depth_normal_column) & !is.na(depth_tumor_column)){
        DP<-as.numeric(f[depth_normal_column]) + as.numeric(f[depth_tumor_column])
      } else if(length(grep("DP=",f))>0){
        DP<-strsplit(strsplit(f[grep("DP=",f)], "DP=")[[1]][2],";")[[1]][1]
      } else if(length(grep("t_alt_count", f))>0 & length(grep("t_ref_count", f))>0){
        alt<-strsplit(strsplit(f[grep("t_alt_count", f)],"t_alt_count=")[[1]][2],";|,|_")[[1]][1]
        ref<-strsplit(strsplit(f[grep("t_ref_count", f)],"t_ref_count=")[[1]][2],";|,|_")[[1]][1]
        if(!is.null(alt)){
          DP<-as.numeric(ref) + as.numeric(alt)
        }
      }

      #TDP:Tumor Depth
      TDP<-0
      if(!is.na(depth_tumor_column)){
        TDP<-as.numeric(f[depth_tumor_column])
      } else if(length(grep("\\|1:",f))>0){
        TDP<-sum(as.numeric(rev(strsplit(strsplit(f[grep("\\|1:",f)], "\\|1:")[[1]][2],":")[[1]])[-1]))
      }else if(!is.null(alt) & !is.null(ref)){
        TDP<-as.numeric(alt)
      }

      #Mutation Start/End Position
      m_start<-as.numeric(f[mutation_start_column])
      m_end<-as.numeric(f[mutation_end_column])

      #Ref./Alt on Mutation Position
      m_ref<-f[mutation_ref_column]
      m_alt<-f[mutation_alt_bnd_column]
      m_alt_chr<-strsplit(strsplit(m_alt, ":")[[1]][1], "\\]|\\[")[[1]][2]
      m_alt_pos<-strsplit(strsplit(m_alt, ":")[[1]][2], "\\]|\\[")[[1]][1]
      m_alt_dir<-ifelse(nchar(strsplit(m_alt, "\\[|\\]")[[1]][1]) > 0, "head", "tail")
      m_alt_tra<-ifelse(length(grep("\\[", m_alt))==1, "forward", "reverse")
      m_alt_cont<-ifelse(m_alt_dir == "head",
                         strsplit(m_alt, "\\[|\\]")[[1]][1],
                         rev(strsplit(m_alt, "\\[|\\]")[[1]])[1])
      m_alt_cont<-ifelse(m_alt_dir == "head",
                         substr(m_alt_cont, 2, nchar(m_alt_cont)),
                         substr(m_alt_cont, 1, nchar(m_alt_cont) - 1))
      m_alt_cont<-ifelse(m_alt_dir == "head",
                         tolower(m_alt_cont),
                         paste(rev(trans_to[match(strsplit(tolower(m_alt_cont), "")[[1]], trans_from)]), collapse=""))

      #When Including MultipleIDs
      #For example, f[10]=SAMD11:NM_152486:exon9:c.C880T:p.Q294X...
      if(is.na(nm_id_column) | nm_id_column == 0){
        g_name<-f[gene_symbol_column]
        nm_ids<-list_nm_cut[which(!is.na(match(list_nm_gene, g_name)))]
      } else {
        nm_ids<-f[nm_id_column]
      }
      nm_ids<-nm_ids[grep("NM_", nm_ids)]
      #Calculate All NM_IDs in Each Mutation
      for(nm_id in nm_ids){
        #For example, nm_ids[[1]]="SAMD11"    "NM_152486" "exon9"     "c.C880T"   "p.Q294X"...
        #Obtain refFLAT Data
        s_variants<-match(nm_id, list_nm_cut)
        if(is.na(s_variants)) {
          print("No MATCH!!")
          next
        }
        #Calculate Sets for NM_ID, because NM_id:ExonRegion is not unique!!
        for(v in s_variants){
          for(mate_whether_exon_intron in c("exon", "intron")){
            nm_sep<-strsplit(list_nm[v], "\t")[[1]]
            #Skip Such As "ch5_hap"
            if(nchar(nm_sep[3]) > 5) {
              print("Invarid NM_ID")
              next
            }
            strand<-nm_sep[4]

            trans_start<-as.numeric(nm_sep[7])
            trans_end<-as.numeric(nm_sep[8])
            if(trans_start==trans_end) {
              print("Translation Position Violation")
              next
            }
            exon_start<-as.numeric(strsplit(nm_sep[10], ",")[[1]])
            exon_end<-as.numeric(strsplit(nm_sep[11], ",")[[1]])

            #Obtain DNA sequence of Transcriptome
            #DNAseq is Unique
            dna<-list_fl_dna[match(nm_id, list_fl_NMID)]
            if(is.na(dna)){
              print(paste(nm_id, "was not found in refMrn."))
              next
            }
            if(nchar(dna) != sum(exon_end - exon_start)){
              dif<-sum(exon_end - exon_start) - nchar(dna)
              if(dif > 200 | dif < -400) {
                print("Splicing Variants!!")
                next
              }
            }

            #Get Relative Mutation Position
            if(whether_exon_intron == "exon" & length(which((exon_start - m_start) * (exon_end - m_start) <= 0)) == 0){
              print("Annotation is Exonic but Position is Intronic")
              next
            } else if(whether_exon_intron == "intron" & length(which((exon_start - m_start) * (exon_end - m_start) <= 0)) != 0){
              print("Annotation is Intronic but Position is Exonic")
              next
            } else if(m_start < exon_start[1] | exon_end[length(exon_end)] < m_start){
              print("Annotation is Intronic but Position is not Intron")
              next
            }

            if(whether_exon_intron == "intron") {
              #Mutation Position May Be 1-based Position
              if(strand =="+"){
                point <- exon_end < m_start
                m_point_master <- sum((exon_end - exon_start)[point])
                m_point_slave  <- sum((exon_end - exon_start)[point]) + 1
              } else {
                point<- m_start < exon_start
                m_point_master <-sum((exon_end - exon_start)[point])
                m_point_slave  <-sum((exon_end - exon_start)[point]) + 1
              }
            } else if(whether_exon_intron == "exon" & mate_whether_exon_intron == "exon"){
              if(strand=="+"){
                point<-(exon_end > m_start)
                plus<-0
                if(length(which(point))>0){
                  if(m_start - exon_start[which(point)[1]] > 0) {
                    plus<-m_start - exon_start[which(point)[1]]
                  }
                }
                m_point_master <- sum((exon_end - exon_start)[!point]) + plus
                m_point_slave  <- sum((exon_end - exon_start)[!point]) + plus + 1
              } else {
                point<-(exon_start > m_start)
                plus<-0
                if(length(which(!point))>0){
                  if((exon_end[rev(which(!point))[1]] - m_start) + 1 > 0){
                    plus<-(exon_end[rev(which(!point))[1]] - m_start) + 1
                  }
                }
                m_point_master <- sum((exon_end - exon_start)[point]) + plus
                m_point_slave  <- sum((exon_end - exon_start)[point]) + plus + 1
              }
            } else {
              #Mutation Position May Be 1-based Position
              if(strand =="+"){
                point <- exon_end < m_start
                m_point_master <- sum((exon_end - exon_start)[point])
                point <- exon_start < m_start
                m_point_slave  <- sum((exon_end - exon_start)[point])
              } else {
                point<- m_start < exon_start
                m_point_master <-sum((exon_end - exon_start)[point])
                point<- m_start < exon_end
                m_point_slave  <-sum((exon_end - exon_start)[point])
              }
            }

            #Get Relative Translation-Start Position (Correct 0-start to 1-start)
            if(strand== "+"){
              point<-(exon_end > trans_start)
              ts_point<-sum((exon_end - exon_start)[!point]) +
                (trans_start - exon_start[which(point)[1]]) + 1
            }else{
              point<-(exon_start > trans_end)
              ts_point<-sum((exon_end - exon_start)[point]) +
                (exon_end[rev(which(!point))[1]] - trans_end) + 1
            }

            #Check Start Codon
            d<-0
            if(substr(dna, ts_point, ts_point + 2)!="atg"){
              flag<-FALSE
              for(d in -2:2){
                if(substr(dna, ts_point + d, ts_point + 2 + d)=="atg"){
                  flag<-TRUE
                  break
                }
              }
              if(flag){
                if(d < 0){
                  dna<-sub(" ", "", paste(paste(rep("x", -d),collapse=""), dna, collapse=""))
                }else{
                  dna<-substr(dna, d + 1, nchar(dna))
                }
              }else{
                print("Non-Start!!")
                next
              }
            }

            #Get Relative Translation-End Position
            if(strand=="+"){
              point<-(exon_end >= trans_end)
              te_point<-sum((exon_end - exon_start)[!point]) +
                (trans_end - exon_start[which(point)[1]])
            }else{
              point<-(exon_start > trans_start)
              te_point<-sum((exon_end - exon_start)[point]) +
                (exon_end[rev(which(!point))[1]] - trans_start)
            }

            #Check Stop Codon
            e<-0
            if(amino[match(substr(dna, te_point-2, te_point), codon)]!="X"){
              dna_trans<-substr(dna, ts_point, te_point)
              flag_2<-FALSE
              dif<-nchar(dna_trans)%%3
              for(e in (seq(from=-3,to=30,3) - dif)){
                if(nchar(dna) < te_point) break
                if(amino[match(substr(dna, te_point - 2 + e, te_point + e), codon)]=="X"){
                  te_point <- te_point + e
                  flag_2<-TRUE
                  break
                }
              }
              if(!flag_2){
                print("Non-Stop!!")
                next
              }
            }

            #Make Peptide
            dna_full<-dna
            dna_trans<-substr(dna, ts_point, te_point)

            #Make Normal Peptide
            peptide_normal<-NULL
            dna<-dna_trans
            while(nchar(dna_trans)>=3){
              peptide_normal<-c(peptide_normal, amino[match(substr(dna_trans, 1, 3), codon)])
              dna_trans<-substr(dna_trans, 4, nchar(dna_trans))
            }

            if(nchar(dna_trans)%%3!=0) {
              print("Peptide_Length_Miss!!")
              next
            }

            if(strand == "+" & m_alt_dir == "head") {
              master_slave <- "master"
              if(m_alt_tra == "forward") {
                mate_strand <- "+"
              } else {
                mate_strand <- "-"
              }
            }
            if(strand == "+" & m_alt_dir == "tail") {
              master_slave <- "slave"
              if(m_alt_tra == "forward") {
                mate_strand <- "+"
              } else {
                mate_strand <- "-"
              }
            }
            if(strand == "-" & m_alt_dir == "head") {
              master_slave <- "slave"
              if(m_alt_tra == "forward") {
                mate_strand <- "-"
              } else {
                mate_strand <- "+"
              }
            }
            if(strand == "-" & m_alt_dir == "tail") {
              master_slave <- "master"
              if(m_alt_tra == "forward") {
                mate_strand <- "+"
              } else {
                mate_strand <- "-"
              }
            }
            m_point_master_2 <- m_point_master - (ts_point) + 1
            m_point_slave_2  <- m_point_slave  - (ts_point) + 1

            dna_master <- substr(dna, 1, m_point_master_2)
            if(master_slave == "master" &
               whether_exon_intron == "exon" &
               mate_whether_exon_intron == "exon" &
               nchar(m_alt_cont) > 0){
              paste(dna_master, m_alt_cont, sep = "")
            }
            dna_slave  <- substr(dna_full, m_point_slave, nchar(dna_full))
            frame <- nchar(substr(dna, 1, m_point_slave_2 - 1))
            if(frame == 0) frame <- m_point_slave - ts_point

            set<-rbind(set, c(g_name,
                              nm_id,
                              strand,
                              m_ref,
                              m_alt,
                              m_alt_dir,
                              m_alt_tra,
                              exon_start[1],
                              rev(exon_end)[1],
                              m_start,
                              dna_master,
                              dna_slave,
                              paste(peptide_normal, collapse=""),
                              chr,
                              whether_exon_intron,
                              mate_whether_exon_intron,
                              master_slave,
                              frame))
          }
        }
      }
      sets<-c(sets, list(set))
    }

    if(length(sets) < 2) next

    count_main <- 0
    for(main_set in sets){
      count_main <- count_main + 1
      if(is.null(main_set)) next

      count_sub <- 0
      for(sub_set in sets){
        count_sub <- count_sub + 1
        if(is.null(sub_set)) next

        if(count_main == count_sub) next

        for(row1 in 1:nrow(main_set)){
          d1<-main_set[row1, ]
          if(d1[17] =="slave" | nchar(d1[11])==0) next

          for(row2 in 1:nrow(sub_set)){
            d2<-sub_set[row2, ]
            if(d2[17] =="master" | nchar(d2[12])==0) next

            if(d1[16] != d2[15] | d2[16] != d1[16]) next
            if(d1[1] == d2[1] & d1[2] != d2[2]){
              next
            }

            flag<-FALSE
            dna_trans<-""
            if( (d1[3] == "+" & d1[6] == "head" & d1[7] == "forward" & d2[3] == "+") |
                (d1[3] == "+" & d1[6] == "head" & d1[7] == "reverse" & d2[3] == "-") |
                (d1[3] == "-" & d1[6] == "tail" & d1[7] == "forward" & d2[3] == "+") |
                (d1[3] == "-" & d1[6] == "tail" & d1[7] == "reverse" & d2[3] == "-") ){
              dna_trans<-paste(d1[11], d2[12], collapse = "", sep="")
              flag<-TRUE
            }

            if(flag){
              peptide<-NULL
              dna_mut<-dna_trans
              while(nchar(dna_trans) >= 3 & length(grep("X", peptide)) == 0){
                peptide<-c(peptide, amino[match(substr(dna_trans, 1, 3), codon)])
                dna_trans<-substr(dna_trans, 4, nchar(dna_trans))
              }

              #Calculate Peptide Start
              peptide_normal<-strsplit(d1[13], "")[[1]]
              min_len<-min(length(peptide), length(peptide_normal))
              peptide_start<-which(peptide[1:min_len] != peptide_normal[1:min_len])[1] - max_peptide_length + 1
              if(is.na(peptide_start)) break
              if(peptide_start < 1) peptide_start<-1

              #Calculate Peptide End
              if(as.numeric(d2[18]) >= 0){
                frame<-ifelse(nchar(d1[11])%%3 == (as.numeric(d2[18]) %% 3), "In", "Out")
              } else {
                frame<-ifelse((nchar(d1[11]) + abs(as.numeric(d2[18]))) %% 3 == 0, "In", "Out")
              }
              if(frame == "In"){
                peptide_normal<-strsplit(d2[13], "")[[1]]
                min_len<-min(length(peptide), length(peptide_normal))
                peptide_end<-which(rev(peptide)[1:min_len] != rev(peptide_normal)[1:min_len])[1]
                if(is.na(peptide_end)) peptide_end <- min_len
                peptide_end <- length(peptide) - peptide_end + max_peptide_length + 10
                if(peptide_end > length(peptide)) peptide_end = length(peptide)
              } else {
                peptide_end = length(peptide)
              }

              peptide <- peptide[peptide_start:peptide_end]
              peptide_normal <- paste(d1[13], "X", d2[13], sep="", collapse="")
              #Save Peptide
              if(length(grep(paste(peptide, sep="", collapse=""), refFasta[,15])) > 0) next
              X<-grep("X", peptide)
              if(length(X) > 0 & IgnoreShortPeptides){if(X < 8) next}
              if(max_peptide_length >= 15 & length(X) > 0 & IgnoreShortPeptides){if(X < 15) next}
              refFasta<-rbind(refFasta, c(paste(random, d1[1], sep="_"),
                                          d1[14],
                                          paste(d1[2], d2[2], sep="_"),
                                          paste(frame, d1[1], d1[15], d2[1], d2[15], sep="_"),
                                          d1[4],
                                          d1[5],
                                          0,
                                          0,
                                          d1[8],
                                          d1[9],
                                          d1[10],
                                          0,
                                          0,
                                          peptide_normal,
                                          paste(peptide, sep="", collapse=""),
                                          0,
                                          dna_mut))

              #Remove X and Save Fasta
              if(!is.na(match("X",peptide)))
                peptide<-peptide[1:(match("X",peptide) - 1)]
              fasta<-c(fasta, sub("_","",paste(">", random, d1[1], sep="_")))
              fasta<-c(fasta, paste(peptide, collapse=""))
              random<-random + 1
              print(paste("OK", uID))
            }
          }
        }
      }
    }
  }

  #Integrate The Same Peptide
  if(is.null(refFasta)) q("no")
  i<-1
  if(!is.null(nrow(refFasta))){
    while(i<=nrow(refFasta)){
      #if gene symbol, mutation position, mutated peptide, normal peptide
      hit<-which((refFasta[i, 11]==refFasta[,11]) &
                   (refFasta[i, 14]==refFasta[,14]) &
                   (refFasta[i, 15]==refFasta[,15]))
      if(length(hit)==1){
        i<-i+1
        next
      }
      #collapse NM_ID, Info, Exon-start, Exon-end
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
              paste(export_dir, "/", rev(strsplit(input_file, "/")[[1]])[1], ".", job_id, ".", "peptide", ".", "fasta", sep=""),
              row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")
  write.table(refFasta,
              paste(export_dir, "/", rev(strsplit(input_file, "/")[[1]])[1], ".", job_id, ".", "peptide", ".", "txt", sep=""),
              row.names=seq(1:nrow(refFasta)), col.names=FALSE, quote=FALSE, sep="\t")
}
