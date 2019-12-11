GenerateMutatedSeq<-function(input_file,
                             hmdir,
                             job_id,
                             refflat_file,
                             refmrna_file,
                             max_peptide_length,
                             chr_column,
                             mutation_start_column,
                             mutation_end_column,
                             mutation_ref_column,
                             mutation_alt_column,
                             nm_id_column,
                             depth_normal_column,
                             depth_tumor_column,
                             ambiguous_between_exon,
                             ambiguous_codon,
                             export_dir,
                             IgnoreShortPeptides,
                             SNPs = NA,
                             multiple_variants = TRUE,
                             apply_annotation = FALSE){

  #READ Data
  data <- read_data(input_file)
  data <- data[grep("\texonic\t", apply(data, 1, function(x) paste(x, collapse = "\t"))), ]
  data <- data[grep("\tnonsynonymous", apply(data, 1, function(x) paste(x, collapse = "\t"))), ]

  if(nrow(data) < 1 | is.null(data)) return(NULL)

  #READ refFlat
  list_nm <- read_refFlat(refflat_file)
  list_nm_gene <- list_nm[, 1]
  list_nm_cut <- list_nm[, 2]

  #READ SNPs Data if available
  if(!is.na(SNPs)) SNPs_vcf <- read_data(SNPs)

  #Get RNA-Code Data
  list_mra <- read_refmrn(refmrna_file)
  start_ <- grep(">", list_mra)
  end_ <- c(start_[-1] - 1, length(list_mra))
  list_fl_NMID <- gsub(">", "", sapply(list_mra[start_], function(x) strsplit(x, " ")[[1]][1]))
  list_fl_dna <- sapply(1:length(start_), function(x) paste(list_mra[(start_[x] + 1):end_[x]], collapse = ""))

  trans_from <- c("a", "t", "g", "c")
  trans_to <- c("t", "a", "c", "g")

  fasta <- NULL
  fasta_wt <- NULL
  refFasta <- NULL
  id <- 0
  for(i in 1:nrow(data)){
    print(paste("Start Analysis: Mutation", i))

    #Extract i-th Data
    f <- as.character(data[i, ])

    #Chromosome
    chr<-f[chr_column]

    #MP:Somatic Mutation Probability
    MP<-0
    if(length(grep("MP=", f))>0){
      MP <- as.numeric(strsplit(strsplit(f[grep("MP=",f)], "MP=")[[1]][2],";")[[1]][1])
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
      DP <- as.numeric(f[depth_normal_column]) + as.numeric(f[depth_tumor_column])
    } else if(length(grep("DP=",f))>0){
      DP <- strsplit(strsplit(f[grep("DP=",f)], "DP=")[[1]][2],";")[[1]][1]
    } else if(length(grep("t_alt_count", f))>0 & length(grep("t_ref_count", f))>0){
      alt <- strsplit(strsplit(f[grep("t_alt_count", f)], "t_alt_count=")[[1]][2],";|,|_")[[1]][1]
      ref <- strsplit(strsplit(f[grep("t_ref_count", f)], "t_ref_count=")[[1]][2],";|,|_")[[1]][1]
      if(!is.null(alt)){
        DP <- as.numeric(ref) + as.numeric(alt)
      }
    }

    #TDP:Tumor Depth
    TDP <- 0
    if(!is.na(depth_tumor_column)){
      TDP <- as.numeric(f[depth_tumor_column])
    } else if(length(grep("\\|1:",f))>0){
      TDP <- sum(as.numeric(rev(strsplit(strsplit(f[grep("\\|1:",f)], "\\|1:")[[1]][2],":")[[1]])[-1]))
    }else if(!is.null(alt) & !is.null(ref)){
      TDP <- as.numeric(alt)
    }

    #Mutation Start/End Position
    m_start <- as.numeric(f[mutation_start_column])
    m_end <- as.numeric(f[mutation_end_column])

    #Ref./Alt on Mutation Position
    m_ref <- f[mutation_ref_column]
    m_alt <- f[mutation_alt_column]

    #When Including MultipleIDs
    #For example, f[nm_id_column]=SAMD11:NM_152486:exon9:c.C880T:p.Q294X...
    nm_ids <- strsplit(f[nm_id_column], ":|,|;")
    hit <- as.numeric(sapply(nm_ids, function(x) grep("NM_", x)))

    #Calculate All NM_IDs in Each Mutation
    Pass <- FALSE
    for(h in hit){
      #For example, nm_ids[[1]]="SAMD11"    "NM_152486" "exon9"     "c.C880T"   "p.Q294X"...
      g_name <- nm_ids[[1]][h - 1]
      nm_id <- nm_ids[[1]][h]
      ans <- strsplit(strsplit(nm_ids[[1]][h+3], "\\.")[[1]][2], "[0-9]")[[1]]
      ans_from <- ans[1]
      ans_to <- rev(ans)[1]
      ans_acid <- strsplit(strsplit(nm_ids[[1]][h+2], "\\.")[[1]][2], "[0-9]")[[1]]
      ans_acid_from <- ans_acid[1]
      ans_acid_to <- rev(ans_acid)[1]

      #Obtain refFLAT Data
      s_variants<-match(nm_id, list_nm_cut)
      if(is.na(s_variants)) {
        print(paste("NM_ID NOT Macth, Skip:", nm_id))
        next
      }

      #Calculate Sets for NM_ID, because NM_id:ExonRegion is not unique!!
      for(v in s_variants){
        final_s_variants <- match(v, s_variants) / length(s_variants) == 1
        nm_sep <- sapply(list_nm[v, ], as.character)

        #Skip Such As "ch5_hap"
        if(nchar(nm_sep[3]) > 5) next
        strand <- nm_sep[4]

        #Check Ref/Alt
        if(strand=="+" & m_ref == ans_acid_from){
          print(paste(nm_sep[2], "- Ref Matched to and NM_ID Attached Information"))
        } else if(strand=="-" &
                 tolower(m_ref)==trans_to[match(tolower(ans_acid_from), trans_from)]){
          print(paste(nm_sep[2], "- Ref Matched to and NM_ID Attached Information"))
        } else {
          print(paste("Ref and NM_ID Attached Information Do not Match, Skip:", nm_id))
          next
        }

        #Get Translation Start/End, Exon Start/End
        trans_start<-as.numeric(nm_sep[7])
        trans_end<-as.numeric(nm_sep[8])
        exon_start<-as.numeric(strsplit(nm_sep[10], ",")[[1]])
        exon_end<-as.numeric(strsplit(nm_sep[11], ",")[[1]])

        #Check Whether Mutation is Among Exonic Region
        #0-base(exon_start), 1-based(exon_end, m_start)
        if(length(which(exon_start < m_start & m_start <= exon_end))!=1){
          print(paste("The Mutation is not between Exon Region, Skip", nm_id))
          next
        }

        #Obtain DNA sequence of Transcriptome
        dna <- list_fl_dna[match(nm_id, list_fl_NMID)]

        #Check DNA
        if(check_dna_validity(dna, nm_id, exon_end, exon_start, ambiguous_between_exon, final_s_variants, Pass)) next

        #Get Relative Mutation Position
        m_point <- get_relative_mutation_position(strand, exon_end, m_start, exon_start)

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
        for(k in unique(0, d, e)){
          dna_trans <- substr(dna, ts_point, te_point)
          if(!is.na(SNPs)) dna_trans <- apply_multiple_snps(SNPs_vcf, exon_start, mutation_start_column, exon_end, chr, strand, dna_trans, trans_to, trans_from)
          dna_trans_normal <- dna_trans
          m_point_2 <- m_point - ts_point + 1 - k

          #Mutation Position is not between Translational Region
          if(m_point_2 < 0) {
            print("Mutation Position is not between Translational Region")
            next
          }

          #Translation Region is not Valid
          if(nchar(dna_trans)%%3!=0) {
            print("Translation Region is not Valid.")
            next
          }

          #Make Normal Peptide
          peptide_normal <- make_normal_peptide(dna_trans, amino, codon, k, e)
          if(is.null(peptide_normal)) next
          target_amino_before <- peptide_normal[ceiling(m_point_2 / 3.0)]

          if(match("X", peptide_normal) < length(peptide_normal)){
            print("Invalid peptide was generated.")
            next
          }







          #Apply Multiple SNVs
          dna_trans_mut <- apply_multiple_snvs(data, multiple_variants, i, exon_start, mutation_start_column, exon_end, chr, strand, dna_trans, trans_to, trans_from)

          #Make Mutated-DNA
          dna_trans_mut <- make_mutated_dna(strand, dna_trans_mut, m_point_2, m_ref, m_alt, trans_to, trans_from)
          if(is.null(dna_trans_mut)) next

          #Make Mutated-Peptide
          peptide <- make_mutated_peptide(dna_trans_mut, amino, codon)
          if(is.null(peptide)) next
          target_amino_after <- peptide[ceiling(m_point_2 / 3.0)]

          #VCF Description of Normal Amino Acid is not What Generated
          if(target_amino_before==target_amino_after | target_amino_after == "X"){
            next
          }

          #VCF Description of Mutated Amino Acid is not What Generated
          if(target_amino_after!=ans_to | target_amino_before!=ans_from){
            next
          }

          #Generate Mutated and Normal Peptide
          frac <- generate_fraction(m_point_2, max_peptide_length, peptide)
          if(is.null(frac)) next
          peptide <- peptide[frac]
          peptide_normal <- peptide_normal[frac]

          #Save Peptide
          if(length(peptide) < 8 & IgnoreShortPeptides) break
          if(max_peptide_length >= 15 & IgnoreShortPeptides & length(peptide) < 8) next
          refFasta<-rbind(refFasta,
                         c(paste(id, gsub("\"","", g_name), sep="_"),
                           chr,
                           nm_ids[[1]][h],
                           nm_ids[[1]][h+2],
                           m_ref,
                           m_alt,
                           round(as.numeric(MP),5),
                           ifelse(is.character(GP), GP, round(GP,5)),
                           exon_start[1],
                           rev(exon_end)[1],
                           m_start,
                           DP,
                           TDP,
                           paste(peptide_normal, collapse=""),
                           paste(peptide, collapse=""),
                           dna_trans_normal,
                           dna_trans_mut))

         #Remove X and Save Fasta in Mutated Peptide
         if(!is.na(match("X", peptide))){
           peptide<-peptide[1:(match("X", peptide) - 1)]
         }
         fasta <- c(fasta, sub("_","", paste(">", id, gsub("\"","", g_name), sep="_")))
         fasta <- c(fasta, paste(peptide, collapse=""))

         #Remove X and Save Fasta in Normal Peptide
         if(!is.na(match("X", peptide_normal))){
           peptide_normal<-peptide_normal[1:(match("X", peptide_normal) - 1)]
         }
         fasta_wt<-c(fasta_wt, sub("_", "", paste(">", id, gsub("\"","", g_name), sep="_")))
         fasta_wt<-c(fasta_wt, paste(peptide_normal, collapse=""))

         id <- id + 1
         print("Peptide Successfully Generated!!")
         stop_loop <- TRUE
         Pass <- TRUE
         break
        }
        if(stop_loop) break
      }
    }
    #Notification
    if(!Pass){
      print("refFlat and refMrna data do not Match to vcf Description")
      print(f[1:7])
      print("Skip This Mutation")
    }
  }

  #Integrate The Same Peptide
  if(is.null(refFasta)) return(NULL)

  if(!is.null(nrow(refFasta))){
    integrated_results <- integrate_same_peptide(refFasta, fasta, fasta_wt)
    refFasta <- integrated_results[[1]]
    fasta <- integrated_results[[2]]
    fasta_wt <- integrated_results[[3]]
  }

  write.table(fasta,
              paste(export_dir, "/", rev(strsplit(input_file, "/")[[1]])[1], ".", job_id, ".", "peptide", ".", "fasta", sep=""),
              row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")
  write.table(fasta_wt,
              paste(export_dir, "/", rev(strsplit(input_file, "/")[[1]])[1], ".", job_id, ".", "wtpeptide", ".", "fasta", sep=""),
              row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")
  write.table(refFasta,
              paste(export_dir, "/", rev(strsplit(input_file, "/")[[1]])[1], ".", job_id, ".", "peptide", ".", "txt", sep=""),
              row.names=seq(1:nrow(refFasta)), col.names=FALSE, quote=FALSE, sep="\t")
}
