GenerateMutatedSeq<-function(input_file, hmdir = getwd(), job_ID, 
                             refFlat_file = paste(hmdir,"/../data/refFlat.txt",sep=""), 
                             refMrna_1 = paste(hmdir,"/../data/refMrna.cut1.fa",sep=""), 
                             refMrna_3 = paste(hmdir,"/../data/refMrna.cut3.fa",sep=""),
                             max_peptide_length = 13, Chr_Column = 1, Mutation_Start_Column = 2, 
                             Mutation_End_Column = 3, Mutation_Ref_Column = 4, Mutation_Alt_Column = 5, 
                             NM_ID_Column = 10, Depth_Normal_Column = NA, Depth_Tumor_Column = NA,
                             ambiguous_between_exon = 0, ambiguous_codon = 0){

  #READ Data
  index<-strsplit(scan(input_file, "character", sep="\n", nlines=1), "\t")[[1]]
  data<-scan(input_file, "character", sep="\n", skip=1)
  data<-data[grep("\texonic\t", data)]
  data<-data[grep("\tnonsynonymous", data)]
  data<-gsub("\"", "",data)
  if(length(data)<1) reutn(NULL)

  #READ refFlat
  list_nm<-scan(refFlat_file, "character", sep="\n")
  list_nm_cut<-sapply(list_nm, function(x) strsplit(x, "\t")[[1]][2])
  
  #Get RNA-Code Data
  list_fl_NMID<-scan(refMrna_1, "character", sep="\n")
  list_fl_dna<- scan(refMrna_3, "character", sep="\n")

  trans_from<-c("a", "t", "g", "c")
  trans_to<-c("t", "a", "c", "g")

  pep_len<-max_peptide_length
  fasta<-NULL
  fasta_norm<-NULL
  refFasta<-NULL
  random<-0
  for(i in 1:length(data)){
   print(paste("Start Analysis: Mutation", i))
    
   #Extract i-th Data
   f<-strsplit(data[i], "\t")[[1]]
   
   #Chromosome
   chr<-f[Chr_Column]

   #MP:Somatic Mutation Probability
   MP<-0
   if(length(grep("MP=",f))>0){
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
   if(!is.na(Depth_Normal_Column) & !is.na(Depth_Tumor_Column)){
    DP<-as.numeric(f[Depth_Normal_Column]) + as.numeric(f[Depth_Tumor_Column])
   } else if(length(grep("DP=",f))>0){
    DP<-strsplit(strsplit(f[grep("DP=",f)], "DP=")[[1]][2],";")[[1]][1]
   } else if(length(grep("t_alt_count", f))>0){
    alt<-strsplit(strsplit(f[grep("t_alt_count", f)],"t_alt_count=")[[1]][2],";|,|_")[[1]][1]
    ref<-strsplit(strsplit(f[grep("t_ref_count", f)],"t_ref_count=")[[1]][2],";|,|_")[[1]][1]
    if(!is.null(alt)){
     DP<-as.numeric(ref) + as.numeric(alt)
    }
   }

   #TDP:Tumor Depth      
   TDP<-0
   if(!is.na(Depth_Tumor_Column)){
    TDP<-as.numeric(f[Depth_Tumor_Column])
   } else if(length(grep("\\|1:",f))>0){
    TDP<-sum(as.numeric(rev(strsplit(strsplit(f[grep("\\|1:",f)], "\\|1:")[[1]][2],":")[[1]])[-1]))
   }else if(!is.null(alt)){
    TDP<-as.numeric(alt)
   }

   #Mutation Start/End Position
   m_start<-as.numeric(f[Mutation_Start_Column])
   m_end<-as.numeric(f[Mutation_End_Column])
   
   #Ref./Alt on Mutation Position
   m_ref<-f[Mutation_Ref_Column]
   m_alt<-f[Mutation_Alt_Column]
   
   #When Including MultipleIDs
   #For example, f[NM_ID_Column]=SAMD11:NM_152486:exon9:c.C880T:p.Q294X...
   nm_ids<-strsplit(f[NM_ID_Column], ":|,|;")
   hit<-as.numeric(sapply(nm_ids, function(x) grep("NM_", x)))
   
   #Calculate All NM_IDs in Each Mutation
   Pass<-FALSE
   for(h in hit){
      #For example, nm_ids[[1]]="SAMD11"    "NM_152486" "exon9"     "c.C880T"   "p.Q294X"...
      g_name<-nm_ids[[1]][h - 1]
      nm_id<-nm_ids[[1]][h]
      ans<-strsplit(strsplit(nm_ids[[1]][h+3],"\\.")[[1]][2],"[0-9]")[[1]]
      ans_from<-ans[1]
      ans_to<-rev(ans)[1]
      ans_acid<-strsplit(strsplit(nm_ids[[1]][h+2],"\\.")[[1]][2],"[0-9]")[[1]]
      ans_acid_from<-ans_acid[1]
      ans_acid_to<-rev(ans_acid)[1]

      #Obtain refFLAT Data
      s_variants<-match(nm_id, list_nm_cut)
      if(is.na(s_variants)) {
        print(paste("No Macth, Skip", nm_id))
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

        #Check Ref/Alt
        if(strand=="+" & 
           m_ref==ans_acid_from){
         #print("Ref and NM_ID Attached Information Match")
        }else if(strand=="-" & 
                 tolower(m_ref)==trans_to[match(tolower(ans_acid_from), trans_from)]){
         #print("Ref and NM_ID Attached Information Match")
        } else {
          print(paste("Ref and NM_ID Attached Information Do not Match, Skip", nm_id))
          next
        }
         
         #Get Translation Start/End, Exon Start/End
         trans_start<-as.numeric(nm_sep[7])
         trans_end<-as.numeric(nm_sep[8])
         exon_start<-as.numeric(strsplit(nm_sep[10], ",")[[1]])
         exon_end<-as.numeric(strsplit(nm_sep[11], ",")[[1]])
         
        #Check Whether Mutation is Among Exon Region
         if(length(which(exon_start < m_start & m_start <= exon_end))!=1){
           print(paste("The Mutation is not between Exon Region, Skip", nm_id))
          next
         }
      
         #Obtain DNA sequence of Transcriptome
        #DNAseq is Unique
         dna<-list_fl_dna[match(nm_id, list_fl_NMID)]
         if(nchar(dna) != sum(exon_end - exon_start)){
            dif<-sum(exon_end - exon_start) - nchar(dna)
           if(abs(dif) < ambiguous_between_exon) {
             if(Last && !Pass){
               print(paste("cDNA Length does not Match to Exon-Start/End Length, Skip", nm_id))
               print(paste("Ambiguous Length", ambiguous_between_exon))
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
            dna_trans<-substr(dna, ts_point, te_point)
            m_point_2<-m_point - (ts_point) + 1
            if(strand=="+"){
             substr(dna_trans, m_point_2, m_point_2)<-tolower(m_alt)
            }else{
              substr(dna_trans, m_point_2, m_point_2)<-trans_to[match(tolower(m_alt),trans_from)]
            }

            #Make Mutated-Peptide
            peptide<-NULL
           dna_trans_mut<-dna_trans
            while(nchar(dna_trans)>=3){
              peptide<-c(peptide, amino[match(substr(dna_trans, 1, 3), codon)])
             dna_trans<-substr(dna_trans, 4, nchar(dna_trans))
            }
            target_amino_after<-peptide[ceiling(m_point_2/3.0)]
            
            #VCF Description of Normal Amino Acid is not What Generated
            if(target_amino_before==target_amino_after){
             next
            }
            
            #VCF Description of Mutated Amino Acid is not What Generated
            if(target_amino_after!=ans_to | target_amino_before!=ans_from){
             next
            }
            
            #Generate Mutated and Normal Peptide
            peptide_start<-ceiling(m_point_2/3.0) - pep_len
            if(peptide_start<1) peptide_start<-1
            peptide_end<-ceiling(m_point_2/3.0) + pep_len
            if(peptide_end>length(peptide)) peptide_end<-length(peptide)
            peptide<-peptide[peptide_start:peptide_end]

           #Save Peptide
            refFasta<-rbind(refFasta, 
                           c(paste(random, gsub("\"","", g_name), sep="_"), 
                             chr, 
                             nm_ids[[1]][h], 
                             nm_ids[[1]][h+2], 
                             m_ref, 
                             m_alt, 
                             MP, 
                             GP, 
                             exon_start[1], 
                             rev(exon_end)[1], 
                             m_start, 
                             DP, 
                             TDP,
                             paste(peptide_normal[peptide_start:peptide_end], collapse=""),
                            paste(peptide, collapse=""), dna_trans_normal, dna_trans_mut)
                            )

           #Remove X and Save Fasta in Mutated Peptide
            if(!is.na(match("X",peptide))){
               peptide<-peptide[1:(match("X",peptide) - 1)]
            }
            fasta<-c(fasta, sub("_","", paste(">", random, gsub("\"","", g_name), sep="_")))
            fasta<-c(fasta, paste(peptide, collapse=""))

            #Remove X and Save Fasta in Normal Peptide
            if(!is.na(match("X",peptide_normal))){
               peptide_normal<-peptide_normal[1:(match("X",peptide_normal) - 1)]
            }
           fasta_norm<-c(fasta_norm, sub("_","", paste(">", random, gsub("\"","", g_name), sep="_")))
            fasta_norm<-c(fasta_norm, paste(peptide_normal[peptide_start:peptide_end], collapse=""))

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
     #If gene symbol, mutation position, mutated peptide, normal peptide are all the same, Integrate These Peptides
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
     fasta_norm<-fasta_norm[-(c(hit[-1] * 2 - 1, hit[-1] * 2))]
     i<-i+1
   }
  }
  write.table(fasta, paste(input_file,job_ID,"peptide","fasta",sep="."),
       row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")
  write.table(fasta_norm, paste(input_file,job_ID,"normpeptide","fasta",sep="."),
       row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")
  write.table(refFasta, paste(input_file, job_ID,"peptide","txt",sep="."),
       row.names=seq(1:nrow(refFasta)), col.names=FALSE, quote=FALSE, sep="\t")
}
