root<-"Output.Snv.2"
pse<- "Output.Pse.2"
files<-list.files(root)
black<-scan("./../icgc_pan_next/BlackList1.3.txt","character",sep="\n")

post<-c("MEL","MEL","LUC","LUC","sLUC-SCC","sLUC-AD","BOC","BOC","COC","COC","MAC",
"RCC","RCC","RCC","RCC","HNC","HNC","GAC","GAC","ESC","PRC","PRC","HCC","HCC","HCC",
"HCC","PAC","PAC-END","BLC","THC","BTC","OVC","UTC","CEC","BRC","BRC","BRC",
"sBRC-GBM","sBRC-GBM","sBRC-MED","ML","ML","CLL","AML","CMDI")

pre<-c("SKCM","MELA","LUSC","LUAD","LUSC","LUAD","BOCA","SARC","COAD","READ","BRCA",
"RECA","KIRC","KIRP","KICH","ORCA","HNSC","STAD","GACA","ESAD","PRAD","EOPC","LINC",
"LIHC","LICA","LIRI","PACA","PAEN","BLCA","THCA","BTCA","OV","UCEC","CESC","GBM",
"LGG","PBCA","GBM","LGG","PBCA","DLBC","MALY","CLLE","LAML","CMDI")

cinfo<-c("", "Gene ID", "Chr", "NM_ID", "Change", "ref", "alt", "Prob", "Mutation Prob.",
                "Exon Start", "Exon End", "Mutation Position", "Depth", "TumorDepth",
                "Peptide Normal", "Peptide Mutation", "DNA_Normal", "DNA_Mut",
                "TotalRNA", "TumorRNARatio", "TumorRNA", "nA", "nB", "Checker",
                "MutRatio", "MutRatio Min", "MutRatio Max")

remove<-NULL
Result<-NULL
for(file in files){
   if(length(grep(file, black))>0) next

   #Check File Exist
   fs<-list.files(paste(root, file, sep="/"))
   if(length(fs) == 0){
    remove<-c(remove, paste(root, file, sep="/"))
    next
   }
   fp<-list.files(paste(pse,  file, sep="/"))
   if(length(fp) == 0){
    remove<-c(remove, paste(pse, file, sep="/"))
    next
   }

   #Get HLA Types
   hit<-grep("HLAtype.txt", fs)
   f<-fs[hit]
   if(length(f)==0) {
    remove<-c(remove, paste(root, file, sep="/"))
    next
   }
   hla<-strsplit(scan(paste(root, file, f, sep="/"), "character", sep="\n")[3],"-")[[1]][1]

   #Synonym/Nonsynonym in Exon
   hit<-grep("summary.txt", fs)
   f<-fs[hit]
   if(length(f)==0) {
    remove<-c(remove, paste(root, file), sep="/")
    next
   }
   sy<-scan(paste(root, file, f, sep="/"), "character", sep="\n")
   num_sy<-as.numeric(strsplit(gsub("[ ]+", "\t",sy[match("exonic - synonymous", sy)+1]),"\t")[[1]][2])
   num_ny<-as.numeric(strsplit(gsub("[ ]+", "\t",sy[match("exonic - nonsynonymous", sy)+1]),"\t")[[1]][2])

   #Synonym/Nonsynonym in Pse
   hit<-grep("summary.txt", fp)
   f<-fp[hit]
   if(length(f)==0) {
    remove<-c(remove, paste(pse, file, sep="/"))
    next
   }
   sy<-scan(paste(pse, file, f, sep="/"), "character", sep="\n")
   num_sy_pse<-as.numeric(strsplit(gsub("[ ]+", "\t",sy[match("exonic - synonymous", sy)+1]),"\t")[[1]][2])
   num_ny_pse<-as.numeric(strsplit(gsub("[ ]+", "\t",sy[match("exonic - nonsynonymous", sy)+1]),"\t")[[1]][2])

   #Make Matrix between Synonymous and Non-synonymous
   table=matrix(c(num_sy, num_ny, num_sy_pse, num_ny_pse), nrow=2, byrow=T, 
      dimnames=list(c("SNV", "PSE"), c("Syn","Non-Syn")))
   ft<-fisher.test(table, alternative="greater")

   #Get Num. NeoAntigen on Exon
   hit<-grep("count.txt", fs)
   f<-fs[hit]
   if(length(f)==0) {
    remove<-c(remove, paste(root, file, sep="/"))
    next
   }
   nn<-sapply(scan(paste(root, file, f, sep="/"), "character", sep="\n"), function(x) strsplit(x,"\t")[[1]])
   tag<-as.character(nn[1,])
   NumNeo<-as.character(nn[nrow(nn),])
   if(length(NumNeo)<100) {
    remove<-c(remove, paste(root, file, sep="/"))
    next
   }

   #Get Num. NeoAntigen on Pse
   hit<-grep("count.txt", fp)
   f<-fp[hit]
   if(length(f)==0) {
    remove<-c(remove, paste(pse, file, sep="/"))
    next
   }
   nn<-sapply(scan(paste(pse, file, f, sep="/"), "character", sep="\n"), function(x) strsplit(x,"\t")[[1]])
   if(is.null(dim(nn)))next
   tag_pse<-as.character(nn[1,])
   NumNeoPse<-as.character(nn[nrow(nn),])
   if(length(NumNeoPse)<100) {
    remove<-c(remove, paste(pse, file, sep="/"))
    next
   }

   #Make Matrix and Calculate Dis
   tmp<-NULL
   Len<-as.numeric(NumNeo[4])
   Len_list<-grep("Length",tag)

   index<-grep("Len|rank|IC50", tag)
   #index<-index[grep("Len|wt500", tag[index])]
   #index<-index[grep("Len|FPKMUP-1-", tag[index])]
   #index<-index[grep("Len|Mer", tag[index])]
   for(j in index){
    if(!is.na(match(j, Len_list))){
     Len<-as.numeric(NumNeo[j])
     #if(Len==0) break
     next
    }
    #GetControl
    key0<-ifelse(length(grep("FPKMNO", tag[j]))>0, "FPKMNO", ifelse(length(grep("FPKMUP",tag[j]))>0, "FPKMUP", "FPKMDOWN"))
    if(key0=="FPKMDOWN") next
    key1<-ifelse(length(grep("IC50", tag[j]))>0, "IC50", "rank")
    key2<-ifelse(length(grep("wt", tag[j]))>0, "wt", "")
    key3<-ifelse(length(grep("Mer", tag[j]))>0, "Mer", "Sep")
    key4<-ifelse(length(grep("CCFP-0-", tag[j]))>0, "CCFP-0-", "CCFP-0.8-")
    if(key1=="IC50"){
     if(key2=="") ques<-paste("IC50", c("-50-","-100","-200-","-500-"), sep="")
     if(key2=="wt") ques<-paste("IC50_wt500", c("-50-","-100-","-200-","-500-"), sep="")
     key5<-ifelse(length(grep(ques[1], tag[j]))==1, ques[1],
           ifelse(length(grep(ques[2], tag[j]))==1, ques[2],
           ifelse(length(grep(ques[3], tag[j]))==1, ques[3], ques[4])))
    }else{
     if(key2=="") ques<-paste("rank", c("-0.1-","-0.5","-2-"), sep="")
     if(key2=="wt") ques<-paste("rank_wt2", c("-0.1-","-0.5","-2-"), sep="")
     key5<-ifelse(length(grep(ques[1], tag[j]))==1, ques[1],
           ifelse(length(grep(ques[2], tag[j]))==1, ques[2], ques[3]))
    }
    
    #Nonsynonymous/NeoAntigen, SNV/PSE
    control_1<-intersect(intersect(intersect(intersect(grep(key1,tag),grep(key2,tag)),grep(key3,tag)),
		grep(key4,tag)),grep(key5,tag))
    control_1<-control_1[grep("FPKMNO",tag[control_1])]
    if(length(control_1)!=1) {print("WARNING")}
    #if(as.numeric(NumNeoPse[4])==0) next
    Len_ctrl<-as.numeric(NumNeoPse[Len_list[rev(which(Len_list < control_1))[1]]])
    table2=matrix(c(as.numeric(NumNeo[j]), Len,
    		    as.numeric(NumNeoPse[control_1]), Len_ctrl), 
		    nrow=2, byrow=T,
    		    dimnames=list(c("SNV", "PSE"), c("NeoAntigen","Nosyno")))
    ft1<-fisher.test(table2, alternative="less")

    tmp<-rbind(tmp, c(file, post[match(hla, pre)], 
    		      ft$p.value, ft$estimate,
    		      ft1$p.value, ft1$estimate,
                      as.numeric(NumNeo[j]), Len, 
		      as.numeric(NumNeoPse[control_1]), as.numeric(NumNeoPse[4]))
	      )
    rownames(tmp)[nrow(tmp)]<-paste(tag[j],tag[control_1],sep="-")
   }
   #if(is.null(tmp))next
   Result<-rbind(Result, cbind(rownames(tmp), 
   			  t(apply(tmp[,c(7,8,9,10,3,4,5,6)], 1,
           		   function(x) c(file, post[match(hla, pre)], num_sy, num_ny, num_sy_pse, num_ny_pse, x)))))
}
colnames(Result)<-c("CutOff","File","Tumor",
		    "ExonSyno","ExonNonSyno","PseSyno","PseuNonsyno",
		    "ExonNonSynoNeo","ExonValidNonSyno",
		    "PseNonSynoNeo", "PseValidNonSyno",
		    "P-SynoNonsyno","Odds-SynoNonsyno",
		    "P-NonsynoNeoanti-SNV/PSE","Odds-NonsynoNeoanti-SNV/PSE"
		   )

pdf("All.Class2.Count.Pse.pdf")
p.val.list<-NULL
p.val.list.tag<-NULL
for(cutoff in unique(Result[,1])){
 hit<-which(!is.na(match(Result[,1],cutoff)))
 for(u in unique(Result[,3])){
  hit2<-hit[which(!is.na(match(Result[hit,3],u)))]
  p<-grep("P-",colnames(Result))
  for(j in 1:length(p)){
   hit3<-hit2[as.numeric(Result[hit2, match("PseuNonsyno", colnames(Result))])>5 & 
    as.numeric(Result[hit2, match("ExonNonSyno", colnames(Result))]) > 5]
   if(length(hit3) == 0) next
   hit3<-hit2
   p.val<- 1-pchisq(df = 2 * length(hit3), -2 * sum(log(as.numeric(Result[hit3,p[j]]))))
   p.val.list<-c(p.val.list, p.val)
   p.val.list.tag<-c(p.val.list.tag, paste(cutoff, u, colnames(Result)[p[j]]))
   hist(as.numeric(Result[hit3,p[j]]),
        main=paste(Result[hit3[1],1], "\n", Result[hit3[1],3], colnames(Result)[p[j]], p.val, sep="-"),breaks=10)
  }
 }
}
dev.off()
names(p.val.list)<-unique(Result[,1])
write.table(Result, "All.Class2.Count.Pse.txt", row.names=FALSE, col.names=TRUE,quote=FALSE,sep="\t")
write.table(remove, "RmList.txt", row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")
