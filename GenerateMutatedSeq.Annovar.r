ls(all=T)
character(0)
rm(list=ls(all=TRUE))

codon<-c("TTT","TTC","TTA","TTG","TCT","TCC","TCA","TCG","TAT","TAC","TAA",
"TAG","TGT","TGC","TGA","TGG","CTT","CTC","CTA","CTG","CCT","CCC","CCA","CCG",
"CAT","CAC","CAA","CAG","CGT","CGC","CGA","CGG","ATT","ATC","ATA","ATG","ACT",
"ACC","ACA","ACG","AAT","AAC","AAA","AAG","AGT","AGC","AGA","AGG","GTT","GTC",
"GTA","GTG","GCT","GCC","GCA","GCG","GAT","GAC","GAA","GAG","GGT","GGC","GGA","GGG")
codon<-tolower(codon)

amino<-c("F","F","L","L","S","S","S","S","Y","Y","X","X","C","C","X","W",
"L","L","L","L","P","P","P","P","H","H","Q","Q","R","R","R","R","I","I","I","M",
"T","T","T","T","N","N","K","K","S","S","R","R","V","V","V","V","A","A","A","A",
"D","D","E","E","G","G","G","G")

data<-scan(commandArgs(TRUE)[1], "character", sep="\n")
index<-strsplit(scan(commandArgs(TRUE)[1], "character", sep="\n", nlines=1), "\t")[[1]]

hmdir<-commandArgs(TRUE)[2]
list_nm_cut<-scan(paste(hmdir,"/refFlat.cut.txt",sep=""), "character", sep="\n")
list_nm    <-scan(paste(hmdir,"/refFlat.txt",    sep=""), "character", sep="\n")

#Get RNA-Code Data
list_fl_NMID<-scan(paste(hmdir,"/refMrna.merge.cut1.fa",sep=""), "character", sep="\n")
list_fl_dna<- scan(paste(hmdir,"/refMrna.merge.cut3.fa",sep=""), "character", sep="\n")

trans_from<-c("a", "t", "g", "c")
trans_to<-c("t", "a", "c", "g")

fasta<-NULL
fasta_norm<-NULL
refFasta<-NULL
random<-0
for(i in 1:length(data)){
   gc();gc();
   #Extract i-th Data
   #f is a row of mutation data
   f<-strsplit(data[i], "\t")[[1]]
   chr<-f[1]
   #MP:Somatic Mutation Probability, GP:Genotype Probability, DP:Full Depth, TDP:Tumor Depth
   if(length(grep("MP=",f))>0){
    MP<-as.numeric(strsplit(strsplit(f[grep("MP=",f)], "MP=")[[1]][2],";")[[1]][1])
   }else if(length(grep("medianVAF=",f))>0){
    MP<-as.numeric(strsplit(strsplit(f[grep("medianVAF=",f)], "medianVAF=")[[1]][2],";")[[1]][1])
   }else{
    MP<-0
   }

   if(length(grep("GP=",f))>0){
    GP<-strsplit(strsplit(f[grep("GP=",f)], "GP=")[[1]][2],";")[[1]][1]
   }else{GP<-0}
   
   if(length(grep("DP=",f))>0){
    DP<-strsplit(strsplit(f[grep("DP=",f)], "DP=")[[1]][2],";")[[1]][1]
   }else{DP<-0}
   
   if(length(grep("\\|1:",f))>0){
    TDP<-sum(as.numeric(rev(strsplit(strsplit(f[grep("\\|1:",f)], "\\|1:")[[1]][2],":")[[1]])[-1]))
   }else{TDP<-0}

   #Mutation Start/End Position
   m_start<-as.numeric(f[2])
   m_end<-as.numeric(f[3])
   #Ref./Alt on Mutation Position
   m_ref<-f[4]
   m_alt<-f[5]
   #Skip if not Exon
   if(f[9]=="unknown"){next
   }else{print(i)}
   
   #When Including MultipleIDs
   #For example, f[10]=SAMD11:NM_152486:exon9:c.C880T:p.Q294X...
   nm_ids<-strsplit(f[10], ":|,")
   hit<-as.numeric(sapply(nm_ids, function(x) grep("NM_", x)))
   Pass<-FALSE
   #Calculate All NM_IDs in Each Mutation
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
         print("No MUTCH!!")
         next
      }
      #Calculate Sets for NM_ID, because NM_id:ExonRegion is not unique!!
      for(v in s_variants){
         Last<-match(v, s_variants)/length(s_variants)==1
         nm_sep<-strsplit(list_nm[v], "\t")[[1]]
	 #Skip Such As "ch5_hap"
	 if(nchar(nm_sep[3]) > 5) next
         strand<-nm_sep[4]

	 #Check Variant's Ref/Alt
	 if(strand=="+" & m_ref==ans_acid_from){
	 }else if(strand=="-" & 
	    tolower(m_ref)==trans_to[match(tolower(ans_acid_from), trans_from)]){
	 }else{
	    if(Last && !Pass){
	       print("Ref/Alt are not Match to vcf Description")
	       print(list_nm_cut[s_variants])
	       print(strand)
	       print(c(m_ref, m_alt, ans_acid_from, ans_acid_to))
	       print("Skip")
	    }
	    next
	 }
         
         trans_start<-as.numeric(nm_sep[7])
         trans_end<-as.numeric(nm_sep[8])
         exon_start<-as.numeric(strsplit(nm_sep[10], ",")[[1]])
         exon_end<-as.numeric(strsplit(nm_sep[11], ",")[[1]])

	 #Check Whether Mutation is Among Exon Region
         if(length(which(exon_start < m_start & m_start <= exon_end))!=1){
	    if(Last && !Pass){
               print("Not Among Exon!!")
	       print(nm_id)
	       print(m_start)
	       print(exon_start)
	       print(exon_end)
	       print("Skip")
	    }
	    next
         }
      
         #Obtain DNA sequence of Transcriptome
	 #DNAseq is Unique
         dna<-list_fl_dna[match(nm_id, list_fl_NMID)]
         if(nchar(dna) != sum(exon_end - exon_start)){
            dif<-sum(exon_end - exon_start) - nchar(dna)
	    if(dif > 200 | dif < -400) {
	       if(Last && !Pass){
                  print("Splicing Variants!!")
	          print(nm_id)
		  print(paste(nchar(dna), sum(exon_end - exon_start)))
		  print("Skip")
	       }
	       next
	    }
         }
	 
         #Get Relative Mutation Position
         if(strand=="+"){
            point<-(exon_end > m_start)
	    plus<-0
	    if(length(which(point))>0){
	       if(m_start - exon_start[which(point)[1]] > 0) plus<-m_start - exon_start[which(point)[1]]
	    }
            m_point<-sum((exon_end - exon_start)[!point]) + plus
         }else{
            point<-(exon_start > m_start)
	    plus<-0
	    if(length(which(!point))>0){
	       if((exon_end[rev(which(!point))[1]] - m_start) + 1 > 0) 
	    		plus<-(exon_end[rev(which(!point))[1]] - m_start) + 1
	    }
	    m_point<-sum((exon_end - exon_start)[point]) + plus
         }

         #Get Relative Translation-Start Position (Correct 0-start to 1-start)
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
	       if(Last && !Pass){
                  print("Non-Start!!")
		  print(nm_id)
	          print(substr(dna, ts_point-4, ts_point + 5))
                  print("Skip")
	       }
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
 	    for(e in (seq(from=-3,to=30,3)-dif)){
	       if(nchar(dna) < te_point) break
               if(amino[match(substr(dna, te_point - 2 + e, te_point + e), codon)]=="X"){
                  te_point <- te_point + e
		  flag_2<-TRUE
                  break
               }
	    }
            if(!flag_2){
	       if(Last && !Pass){
	          print("Non-Stop!!")
		  print(nm_id)
    	          print(substr(dna, te_point-4, te_point + 6))
                  print("Skip")
	       }
	       next
	    }
         }

         #Check Peptide Length
	 stop_loop<-FALSE
         for(k in unique(0, d, e)){
            dna_trans<-substr(dna, ts_point, te_point)
            m_point_2<-m_point - (ts_point) + 1 - k
	    if(m_point_2 < 0) {
	       next
	    }
            if(nchar(dna_trans)%%3!=0) {
	       if(Last && !Pass){
	          print("Peptide_Length_Miss!!")
                  print("Skip")
	       }
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
            if(strand=="+"){substr(dna_trans, m_point_2, m_point_2)<-tolower(m_alt)
            }else{substr(dna_trans, m_point_2, m_point_2)<-
      	       trans_to[match(tolower(m_alt),trans_from)]
            }

            #Make Mutated-Peptide
            peptide<-NULL
	    dna_trans_mut<-dna_trans
            while(nchar(dna_trans)>=3){
               peptide<-c(peptide, amino[match(substr(dna_trans, 1, 3), codon)])
	       dna_trans<-substr(dna_trans, 4, nchar(dna_trans))
            }
            target_amino_after<-peptide[ceiling(m_point_2/3.0)]

            if(target_amino_before==target_amino_after){
	       if(Last && !Pass){
	          print("Un-Changed!!")
		  print(c(m_start, m_point, m_point_2, m_ref, m_alt))
		  print(nm_sep)
		  print(list_nm_cut[s_variants])
		  print(c("target_amino_before", "ans_from", "target_amino_after", "ans_to"))
		  print(c(target_amino_before, ans_from, target_amino_after, ans_to))
	       }
	       next
            }
            if(target_amino_after!=ans_to | target_amino_before!=ans_from){
 	       if(Last && !Pass){
	          print("Wrong Change!!")
		  print(c(m_start, m_point, m_point_2, m_ref, m_alt))
		  print(nm_sep)
		  print(list_nm_cut[s_variants])
		  print(c(target_amino_before, ans_from, target_amino_after, ans_to))
	       }
	       next
            }
            peptide_start<-ceiling(m_point_2/3.0) - 13
            if(peptide_start<1) peptide_start<-1
            peptide_end<-ceiling(m_point_2/3.0) + 13
            if(peptide_end>length(peptide)) peptide_end<-length(peptide)
            peptide<-peptide[peptide_start:peptide_end]

	    #Save Peptide
            refFasta<-rbind(refFasta, 
	       c(paste(g_name, random, sep="_"), chr, nm_ids[[1]][h], nm_ids[[1]][h+2], m_ref, m_alt, MP, GP, 
	          exon_start[1], rev(exon_end)[1], m_start, DP, TDP,
	          paste(peptide_normal[peptide_start:peptide_end], collapse=""),
		     paste(peptide, collapse=""), dna_trans_normal, dna_trans_mut))

	    #Remove X and Save Fasta
            if(!is.na(match("X",peptide)))
               peptide<-peptide[1:(match("X",peptide) - 1)]
            fasta<-c(fasta, sub(":","",paste(">", g_name, random, nm_id, m_point, 
               m_ref, m_alt, target_amino_before, target_amino_after,sep=":")))
            fasta<-c(fasta, paste(peptide, collapse=""))

	    #Remove X and Save Fasta
            if(!is.na(match("X",peptide_normal)))
               peptide_normal<-peptide_normal[1:(match("X",peptide_normal) - 1)]
	    fasta_norm<-c(fasta_norm, sub(":","",paste(">", g_name, random, nm_id, m_point, 
               m_ref, m_alt, target_amino_before, target_amino_after,sep=":")))
            fasta_norm<-c(fasta_norm, paste(peptide_normal[peptide_start:peptide_end], collapse=""))

	    random<-random + 1
	    print("OK")
	    stop_loop=TRUE
	    Pass=TRUE
            break
         }
	 if(stop_loop) break
      }
   }
}

#Integrate The Same Peptide
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
   fasta_norm<-fasta_norm[-(c(hit[-1] * 2 - 1, hit[-1] * 2))]
   i<-i+1
 }
}

write.table(fasta, paste(commandArgs(TRUE)[1],"peptide","fasta",sep="."),
		   row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")
write.table(fasta_norm, paste(commandArgs(TRUE)[1],"normpeptide","fasta",sep="."),
		   row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")
write.table(refFasta, paste(commandArgs(TRUE)[1],"peptide","txt",sep="."),
		   row.names=seq(1:nrow(refFasta)), col.names=FALSE, quote=FALSE, sep="\t")

if(FALSE){
 rna<-scan("refMrna.fa","character",sep="\n")
 rna_one_row<-NULL
 for(i in 1:length(rna)){
   if(length(grep(">", rna[i]))==1){
     s<-paste(c(strsplit(rna[i]," ")[[1]],""), collapse="\t")
     rna_one_row<-c(rna_one_row, sub(">","",s))
   }else{
     rna_one_row[length(rna_one_row)]<-
	paste(rna_one_row[length(rna_one_row)], rna[i], sep="")
   }
 }
 write.table(rna_one_row, "refMrna.merge.fa", row.names=FALSE, col.names=FALSE, 
 	sep="\t", quote=FALSE)
}
