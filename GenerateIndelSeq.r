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

data<-scan(commandArgs(TRUE)[2], "character", sep="\n", nlines=1000)
data1<-scan(commandArgs(TRUE)[2], "character", sep="\n", skip=rev(grep("\\#",data))[1])
data1<-t(sapply(data1, function(x) strsplit(x, "\t")[[1]][c(1,2,7,4,5)]))
data1<-data1[data1[,3]!="LOWSUPPORT",]
if(is.null(nrow(data1))) data1<-t(data1)
if(nrow(data1)==0){q("no")}
data1<-cbind(apply(data1, 1, function(x) paste("\t",x[1],"\t",x[2],"\t",sep="")), data1[,4], data1[,5])

data2<-scan(commandArgs(TRUE)[1], "character", sep="\n", skip=1)
hit<-sapply(data1[,1], function(x) grep(x, data2)[1])
data2<-data2[hit]
unknown<-grep("unknown", data2)
data1<-data1[-unknown,]
data2<-data2[-unknown]
if(is.vector(data1)) data1<-t(data1)
if(nrow(data1)==0){q("no")}

hmdir<-commandArgs(TRUE)[3]
list_nm_cut<-scan(paste(hmdir,"/refFlat.cut.txt",sep=""), "character", sep="\n")
list_nm    <-scan(paste(hmdir,"/refFlat.txt",    sep=""), "character", sep="\n")

#Get RNA-Code Data
list_fl_NMID<-scan(paste(hmdir,"/refMrna.merge.cut1.fa",sep=""), "character", sep="\n")
list_fl_dna<- scan(paste(hmdir,"/refMrna.merge.cut3.fa",sep=""), "character", sep="\n")

trans_from<-c("a", "t", "g", "c")
trans_to<-c("t", "a", "c", "g")

pep_len<-as.numeric(commandArgs(TRUE)[4])
fasta<-NULL
fasta_norm<-NULL
refFasta<-NULL
random<-0
for(i in 1:length(data2)){
   print(i)
   #Extract i-th Data
   #f is a row of mutation data
   f1<-c(strsplit(data1[i,1], "\t")[[1]], data1[i,c(2,3)])
   f2<-strsplit(data2[i], "\t")[[1]]
   f2<-c(f2, f2[1],f2[2])[-c(1,2)]
   #Fisher Test
   #DP<-0
   #TDP<-0
   #tumor_indel = 0
   #normal_indel = 0
   #pin<-grep("GT",f1)
   #label<-strsplit(f1[pin],":")[[1]]
   #DP<-sum(as.numeric(strsplit(f1[pin+1],":")[[1]][grep("PR|NR",label)]))
   #TDP<-sum(as.numeric(strsplit(f1[pin+2],":")[[1]][grep("PR|NR",label)]))
   #NI<-sum(as.numeric(strsplit(f1[pin+1],":")[[1]][grep("PU|NU",label)]))
   #NI<-ifelse(DP < NI, DP, NI)
   #TI<-sum(as.numeric(strsplit(f1[pin+2],":")[[1]][grep("PU|NU",label)]))
   #TI<-ifelse(TDP < TI, TDP, TI)
   DP<-0
   TDP<-0
   if(length(grep("MP=",f2[4]))>0){
    MP<-as.numeric(strsplit(strsplit(f2[4], "MP=")[[1]][2], ";")[[1]][1])
   }else if(length(grep("VAF=",f2[4]))>0){
    MP<-as.numeric(strsplit(strsplit(f2[4], "VAF=")[[1]][2], ";")[[1]][1])
   }else{
    MP<-0
   }

   if(length(grep("GP=",f2[4]))>0){
    GP<-strsplit(strsplit(f[grep("GP=",f2[4])], "GP=")[[1]][2],";")[[1]][1]
   }else{GP<-0}

   #MP<-#fisher.test(matrix(c(TDP - TI, TI, DP - NI, NI),nrow=2,byrow=TRUE))$p
   #GP<-#paste(c(TDP - TI, TI, DP - NI, NI), collapse="/")
   #if(NI > (TDP + DP) * 0.01) next
   #if(MP > 0.05) next
   gc();gc();

   chr<-f2[1]
   #Mutation Start/End Position
   m_start<-as.numeric(f2[2])
   m_end<-as.numeric(f2[3])
   #Ref./Alt on Mutation Position
   m_ref<-f1[4]
   m_alt<-f1[5]
     
   #When Including MultipleIDs
   Pass<-FALSE
   g_name<-f2[6]
   nm_id<-strsplit(f2[5],"\\.")[[1]][1]

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

      trans_start<-as.numeric(nm_sep[7])
      trans_end<-as.numeric(nm_sep[8])
      if(trans_start==trans_end) next
      if(trans_end < m_start) {print("Not Within Translational Region");next;}
      if(trans_start > m_start) {print("Not Within Translational Region");next;}
      exon_start<-as.numeric(strsplit(nm_sep[10], ",")[[1]])
      exon_end<-as.numeric(strsplit(nm_sep[11], ",")[[1]])

      #Check Whether Mutation is Among Exon Region
      #0-base(exon_start), 1-based(exon_end, m_start)
      if(length(which(exon_start + 2 < m_start & m_start <= exon_end))!=1){
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
      #Mutation Position May Be 1-based Position
      if(strand =="+"){
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
	    if((exon_end[rev(which(!point))[1]] - m_start) + 1 > 0) plus<-(exon_end[rev(which(!point))[1]] - m_start) + 1
	 }
	 m_point<-sum((exon_end - exon_start)[point]) + plus
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
	 dna_trans_normal<-dna_trans
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
	 dna_trans_mut<-dna_trans
         m_point_2<-m_point - (ts_point) + 1
	 if(m_point_2 < 4) next
         if(strand == "+"){
	 　 if(substr(dna_trans,m_point_2, m_point_2+nchar(m_ref)-1)==
		   paste(substring(tolower(m_ref),1:nchar(m_ref),1:nchar(m_ref)), collapse="")
            ){
	       dna_trans<-paste(substr(dna_trans, 1, m_point_2-1),
	             paste(sapply(substring(m_alt, 1:nchar(m_alt),1:nchar(m_alt)), 
		        function(x) trans_to[match(tolower(x),trans_from)]), collapse=""),
	             substr(dna_trans, m_point_2, nchar(dna_trans)), sep="")	    
	    } else {next}
         } else {
	 　 if(paste(substring(substr(dna_trans,m_point_2-nchar(m_ref)+1,m_point_2), 
		   1:nchar(m_ref),1:nchar(m_ref)),collapse="")==
		   paste(sapply(rev(substring(tolower(m_ref),1:nchar(m_ref),1:nchar(m_ref))), 
		   function(x) trans_to[match(tolower(x),trans_from)]),collapse="")
            ){
	       dna_trans<-paste(substr(dna_trans, 1, m_point_2-nchar(m_ref)),
	             paste(sapply(rev(substring(m_alt, 1:nchar(m_alt),1:nchar(m_alt))), 
		        function(x) trans_to[match(tolower(x),trans_from)]), collapse=""),
	             substr(dna_trans, m_point_2+1, nchar(dna_trans)), sep="")
	    } else {next}
         }
	 #if(m_start == 39254335)aaa

         #Make Mutated-Peptide
         peptide<-NULL
         while(nchar(dna_trans)>=3){
	    a<-amino[match(substr(dna_trans, 1, 3), codon)]
            peptide<-c(peptide, a)
	    if(a=="X") break
	    dna_trans<-substr(dna_trans, 4, nchar(dna_trans))
         }
	 min_len<-min(length(peptide),length(peptide_normal))

         peptide_start<-which(peptide[1:min_len]!=peptide_normal[1:min_len])[1] - pep_len
	 if(is.na(peptide_start))break
	 if(peptide_start < 1) peptide_start<-1
	 peptide_end<-which(peptide[1:min_len]!=peptide_normal[1:min_len])[1]
	 peptide_length<-length(peptide)
         peptide<-peptide[peptide_start:length(peptide)]

	 #Save Peptide
	 if(length(peptide)<5) break
	 #if(length(peptide_normal[peptide_start:min(length(peptide_normal),peptide_end)])==1)aaa
         refFasta<-rbind(refFasta,
	    c(paste(random, g_name, sep="_"), chr, f2[1], length(peptide), m_ref, m_alt, round(as.numeric(MP),5), ifelse(is.character(GP), GP, round(GP,5)),
	    	  exon_start[1], rev(exon_end)[1],m_start, DP + TDP, TDP,
	          paste(peptide_normal[peptide_start:min(length(peptide_normal),peptide_end)], collapse=""),
		  paste(peptide, collapse=""), dna_trans_normal, dna_trans_mut))

	 #Remove X and Save Fasta
         if(!is.na(match("X",peptide)))
            peptide<-peptide[1:(match("X",peptide) - 1)]
         fasta<-c(fasta, sub("_","",paste(">", random, g_name, sep="_")))
         fasta<-c(fasta, paste(peptide, collapse=""))
	 random<-random + 1

	 print("OK")
	 stop_loop=TRUE
	 Pass=TRUE
         break
      }
      if(stop_loop) break
   }
}

#Integrate The Same Peptide
if(is.null(refFasta))q("no")
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

write.table(fasta, paste(commandArgs(TRUE)[1],"peptide.indel","fasta",sep="."),
		   row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")
write.table(refFasta, paste(commandArgs(TRUE)[1],"peptide.indel","txt",sep="."),
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
