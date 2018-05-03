GenerateIndelSeq<-function(input_file,
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
                           export_dir){

  #READ Data
  index<-strsplit(scan(input_file, "character", sep="\n", nlines=1), "\t")[[1]]
  data<-scan(input_file, "character", sep="\n", skip=1)
  data<-data[grep("\texonic\t", data)]
  data<-data[grep("insertion|deletion", data)]
  data<-gsub("\"", "",data)
  if(length(data)<1) return(NULL)

  #READ refFlat
  list_nm<-scan(refflat_file, "character", sep="\n")
  list_nm_cut<-sapply(list_nm, function(x) strsplit(x, "\t")[[1]][2])

  #Get RNA-Code Data
  list_mra<-scan(refmrna_file, "character", sep=" ")
  start_<-grep(">", list_mra)
  end_<-c(start_[-1] - 1, length(list_mra))
  list_fl_NMID<-gsub(">", "", list_mra[start_])
  list_fl_dna <-sapply(1:length(start_), function(x) paste(list_mra[(start_[x] + 2):end_[x]], collapse = ""))

  trans_from<-c("a", "t", "g", "c")
  trans_to<-c("t", "a", "c", "g")

  fasta<-NULL

  refFasta<-NULL
  random<-0
  for(i in 1:length(data)){
    print(paste("Start Analysis: Mutation", i))

    #Extract i-th Data
    f<-strsplit(data[i], "\t")[[1]]

    #Chromosome
    chr<-f[chr_column]

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
    m_alt<-f[mutation_alt_column]

    #When Including MultipleIDs
    #For example, f[nm_id_column]=SAMD11:NM_152486:exon9:c.C880T:p.Q294X...
    nm_ids<-strsplit(f[nm_id_column], ":|,|;")
    hit<-as.numeric(sapply(nm_ids, function(x) grep("NM_", x)))

    #Calculate All NM_IDs in Each Mutation
    Pass<-FALSE
    for(h in hit){
      #For example, nm_ids[[1]]="SAMD11"    "NM_152486" "exon9"     "c.C880T"   "p.Q294X"...
      g_name<-nm_ids[[1]][h - 1]
      nm_id<-nm_ids[[1]][h]







      #Obtain refFLAT Data
      s_variants<-match(nm_id, list_nm_cut)
      if(is.na(s_variants)) {
        print(paste("NM_ID NOT Macth, Skip:", nm_id))
        next
      }

      #Calculate Sets for NM_ID, because NM_id:ExonRegion is not unique!!
      for(v in s_variants){
        #Whether Last or Not
        Last<-match(v, s_variants)/length(s_variants)==1
        nm_sep<-strsplit(list_nm[v], "\t")[[1]]

        #Skip Such As "ch5_hap"
        if(nchar(nm_sep[3]) > 5) next
        strand<-nm_sep[4]












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
        #DNAseq is Unique
        dna<-list_fl_dna[match(nm_id, list_fl_NMID)]
        if(nchar(dna) != sum(exon_end - exon_start)){
          dif<-sum(exon_end - exon_start) - nchar(dna)
          if(abs(dif) <= ambiguous_between_exon) {
            if(Last && !Pass){
              print(paste("cDNA Length does not Match to Exon-Start/End Length, Skip", nm_id))
            } else {
              print(paste("Permit Ambiguous Exonic Region:", nm_id))
            }
            next
          }
        }

        #Get Relative Mutation Position
        if(strand=="+"){
          point<-(exon_end > m_start)
          plus<-0
          if(length(which(point))>0){
            if(m_start - exon_start[which(point)[1]] > 0) {
              plus<-m_start - exon_start[which(point)[1]]
            }
          }
          m_point<-sum((exon_end - exon_start)[!point]) + plus
        }else{
          point<-(exon_start > m_start)
          plus<-0
          if(length(which(!point))>0){
            if((exon_end[rev(which(!point))[1]] - m_start) + 1 > 0){
              plus<-(exon_end[rev(which(!point))[1]] - m_start) + 1
            }
          }
          m_point<-sum((exon_end - exon_start)[point]) + plus
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
        d<-0
        if(substr(dna, ts_point, ts_point + 2)!="atg"){
          flag<-FALSE
          for(d in  (-1 * ambiguous_codon):(ambiguous_codon)){
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
            print("Permit Ambiguous Codon Start")
          }else{
            print(paste("Start Position is not ATG, Skip", nm_id))
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
          flag<-FALSE
          dif<-nchar(dna_trans)%%3
          for(e in (seq(from=-1 * ambiguous_codon, to=ambiguous_codon, 3) - dif)){
            if(nchar(dna) < te_point) break
            if(amino[match(substr(dna, te_point - 2 + e, te_point + e), codon)]=="X"){
              te_point <- te_point + e
              flag<-TRUE
              break
            }
          }
          if(!flag){
            print(paste("End Position Amino Acid is not X, Skip", nm_id))
            next
          } else {
            print("Permit Ambiguous Codon Start")
          }
        }

        #Check Peptide Length
        stop_loop<-FALSE
        for(k in unique(0, d, e)){
          dna_trans<-substr(dna, ts_point, te_point)
          m_point_2<-m_point - (ts_point) + 1 - k

          #Mutation Position is not between Translational Region
          if(m_point_2 < 0) {
            next
          }

          #Translation Region is not VAlid
          if(nchar(dna_trans)%%3!=0) {
            next
          }

          #Make Normal Peptide
          peptide_normal<-NULL
          dna_trans_normal<-dna_trans
          while(nchar(dna_trans)>=3){
            peptide_normal<-c(peptide_normal,amino[match(substr(dna_trans, 1, 3),codon)])
            dna_trans<-substr(dna_trans, 4, nchar(dna_trans))
          }
          if(k==e & match("X", peptide_normal) < length(peptide_normal)){
            next
          }
          target_amino_before<-peptide_normal[ceiling(m_point_2/3.0)]

          #Make Mutated-DNA
          dna_trans<-substr(dna, ts_point, nchar(dna))
          m_point_2<-m_point - (ts_point) + 1
          if(m_point_2 < 4) next
	        if(strand == "+"){
	          if(m_ref == "-"){
	            #Insertion
              dna_trans<-paste(substr(dna_trans, 1, m_point_2 - 1),
                               paste(sapply(substring(m_alt, 1:nchar(m_alt),1:nchar(m_alt)),
                                            function(x) trans_to[match(tolower(x),trans_from)]), collapse=""),
                               substr(dna_trans, m_point_2, nchar(dna_trans)), sep="")
            } else {
              if(substr(dna_trans, m_point_2, m_point_2 + nchar(m_ref) - 1) ==
                  paste(substring(tolower(m_ref), 1:nchar(m_ref), 1:nchar(m_ref)), collapse="")){

                dna_trans<-paste(substr(dna_trans, 1, m_point_2 - 1),
                                 gsub("NA", "", paste(sapply(substring(m_alt, 1:nchar(m_alt), 1:nchar(m_alt)),
                                                             function(x) trans_to[match(tolower(x),trans_from)]), collapse="")),
                                 substr(dna_trans, m_point_2 + nchar(m_ref), nchar(dna_trans)), sep="")
              } else {
                next
              }
            }
          } else {
          if(m_ref == "-"){
            #Insertion
            dna_trans<-paste(substr(dna_trans, 1, m_point_2),
                            paste(sapply(rev(substring(m_alt, 1:nchar(m_alt), 1:nchar(m_alt))),
                                         function(x) trans_to[match(tolower(x),trans_from)]), collapse=""),
                            substr(dna_trans, m_point_2 + 1, nchar(dna_trans)), sep="")
           } else {
             if(paste(substring(substr(dna_trans, m_point_2 - nchar(m_ref) + 1, m_point_2), 1:nchar(m_ref), 1:nchar(m_ref)), collapse="") ==
                 paste(sapply(rev(substring(tolower(m_ref), 1:nchar(m_ref), 1:nchar(m_ref))), function(x) trans_to[match(tolower(x),trans_from)]), collapse="")){

               dna_trans<-paste(substr(dna_trans, 1, m_point_2 - nchar(m_ref)),
                          gsub("NA","",paste(sapply(rev(substring(m_alt, 1:nchar(m_alt), 1:nchar(m_alt))),
                                                    function(x) trans_to[match(tolower(x),trans_from)]), collapse="")),
                          substr(dna_trans, m_point_2 + 1, nchar(dna_trans)), sep="")
             } else {
               next
             }
           }
         }
         dna_trans_mut<-dna_trans

         #Make Mutated-Peptide
         peptide<-NULL
         while(nchar(dna_trans) >= 3){
           a<-amino[match(substr(dna_trans, 1, 3), codon)]
           peptide<-c(peptide, a)
           if(a=="X") break
           dna_trans<-substr(dna_trans, 4, nchar(dna_trans))
         }












         #Generate Mutated and Normal Peptide
         min_len<-min(length(peptide), length(peptide_normal))
         peptide_start<-which(peptide[1:min_len] != peptide_normal[1:min_len])[1] - max_peptide_length + 1
         if(is.na(peptide_start)) break
         if(peptide_start < 1) peptide_start<-1
         peptide_end<-which(rev(peptide)[1:min_len] != rev(peptide_normal)[1:min_len])[1]
         if(is.na(peptide_end)) peptide_end <- min_len
         peptide_end <- min_len - peptide_end + max_peptide_length + 10
         if(peptide_end > length(peptide)) peptide_end = length(peptide)
         peptide <- peptide[peptide_start:peptide_end]
         peptide_end <- peptide_end + 10
         if(peptide_end > length(peptide_normal)) peptide_end = length(peptide_normal)
         peptide_normal <- peptide_normal[peptide_start:min(peptide_end, length(peptide_normal))]

         #Save Peptide
         X<-grep("X", peptide)
         if(length(X) > 0){if(X < 8) next}
         if(max_peptide_length >= 15 & length(X) > 0){if(X < 15) next}
         frame <- ifelse(abs(nchar(gsub("-", "", m_alt)) - nchar(gsub("-", "", m_ref))) %% 3 == 0, "In", "Out")
         refFasta<-rbind(refFasta,
                         c(paste(random, gsub("\"","", g_name), sep="_"),
                           chr,
                           nm_ids[[1]][h],
                           paste(frame, nm_ids[[1]][h+2], sep="_", collapse="_"),
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

         #Remove X and Save Fasta
         if(!is.na(match("X", peptide))){
           peptide<-peptide[1:(match("X",peptide) - 1)]
         }
         fasta<-c(fasta, sub("_","", paste(">", random, gsub("\"","", g_name), sep="_")))
         fasta<-c(fasta, paste(peptide, collapse=""))








         random<-random + 1
         print("Peptide Successfully Generated!!")
         stop_loop=TRUE
         Pass=TRUE
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
  if(is.null(refFasta)) {
    return(NULL)
  }
  i<-1
  if(!is.null(nrow(refFasta))){
   while(i<=nrow(refFasta)){
     #If mutation position, mutated peptide, normal peptide are all the same, Integrate These Peptides
     #Note That, If the mutation position and chromosome number are the same, Merge script integrate them in the later process.
     hit<-which((refFasta[i, 11]==refFasta[,11]) &
                (refFasta[i, 14]==refFasta[,14]) &
                (refFasta[i, 15]==refFasta[,15]))
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
              paste(export_dir, "/", input_file, ".", job_id, ".", "peptide", ".", "fasta", sep=""),
              row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")
  write.table(refFasta,
              paste(export_dir, "/", input_file, ".", job_id, ".", "peptide", ".", "txt", sep=""),
              row.names=seq(1:nrow(refFasta)), col.names=FALSE, quote=FALSE, sep="\t")
}
