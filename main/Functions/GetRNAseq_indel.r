GetRNAseq_indel<-function(output_peptide_txt_file, RNAseq_file = NA, output_file_rna_vcf = NA){
  data<-t(sapply(scan(output_peptide_txt_file, "character", sep="\n"), function(x) strsplit(x, "\t")[[1]]))
  if(is.na(RNAseq_file) | !file.exists(RNAseq_file)){
     ratio_matrix<-matrix(nrow=nrow(data), ncol=3, NA)
     write.table(cbind(data, ratio_matrix), commandArgs(TRUE)[1],
                 row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")
     q("no")
  }
  
  ##Get RNA Data
  temp<-scan(RNAseq_file, "character", sep="\n")
  if(length(temp) < 2){
    ratio_matrix<-matrix(nrow=nrow(data), ncol=3, NA)
    write.table(cbind(data, ratio_matrix), commandArgs(TRUE)[1],
                row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")
    q("no")
  }
  rna<-t(sapply(temp, function(x) strsplit(x, "\t")[[1]]))
  rna<-rna[-1,]
  rna_pos<-t(sapply(rna[,2], function(x) strsplit(x, ":|-")[[1]]))
  hit<-match(sapply(data[,2], function(x) strsplit(x,"_")[[1]][2]), rna[,1])
  
  #data[,chr], data[,m_position]
  for(x in which(is.na(hit))){
     hit_pos<-which(data[x,3]==rna_pos[,1])
     abs_hit_pos<-hit[which(as.numeric(data[x,12]) > as.numeric(rna_pos[hit,2]) & 
                              as.numeric(data[x,12]) < as.numeric(rna_pos[hit,3]))]
     if(length(abs_hit_pos)!=0){
        print(abs_hit_pos[1])
        hit[x]<-abs_hit_pos[1]
     }
  }
  
  data<-cbind(data, rna[hit,3])
  colnames(data)<-NULL
  rownames(data)<-NULL
  tail_col<-ncol(data)
  
  ##Get Ratio
  if(is.na(output_file_rna_vcf) | !file.exists(output_file_rna_vcf)){
    ratio<-t(sapply(scan(output_file_rna_vcf, "character", sep="\n"),
                    function(x) strsplit(x, "\t")[[1]]))
    ratio<-t(sapply(scan(output_file_rna_vcf, "character", sep="\n", skip=which(lapply(ratio, length)>2)[1]), 
                    function(x) strsplit(x, "\t")[[1]]))
    if(ncol(ratio)==0){
      ratio_matrix<-matrix(nrow=nrow(data), ncol=3, NA)
    } else {
      size<-sum(as.numeric(sapply(strsplit(ratio[1,3],"M|D|I")[[1]], function(x) rev(strsplit(x,"N|S|H|P")[[1]])[1])))
      ratio_matrix<-NULL
      for(i in 1:nrow(data)){
         hit<-which(!is.na(match(ratio[,1], data[i,3])))
         hit<-hit[which(as.numeric(ratio[hit,2]) + size > data[i,12]
            & as.numeric(ratio[hit,2]) -1 < data[i,12])]
         if(length(hit)==0){
            ratio_matrix<-rbind(ratio_matrix, c("0/0.1","0"))
            next
         }
         total<-length(hit) + length(grep("I|D", ratio[hit,3]))
         mut<-length(grep("I|D", ratio[hit,3]))
         r<-paste(mut, total, sep="/")
         ratio_matrix<-rbind(ratio_matrix, c(r, mut/total * as.numeric(data[i,19])))
      }
    }
  }else {
    ratio_matrix<-matrix(nrow=nrow(data), ncol=2, NA)
  }
  write.table(cbind(data, ratio_matrix), output_peptide_txt_file, 
              row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")
}  
