RNAExpression<-function(rnaexp_file,
                        output_peptide_txt_file,
                        width,
                        samtools_dir,
                        refdna_file,
                        rnabam_file,
                        bcftools_dir,
                        indel){
  #Attach RNAseq Data if Exist, Otherwise set NULL column
  skip=FALSE
  if(!is.na(rnaexp_file) & !is.na(rnabam_file) & !is.na(samtools_dir) & !is.na(bcftools_dir)){
    if(!file.exists(rnaexp_file) | !file.exists(rnabam_file)) {
      print("Indicated RNA File or RNA bam file do not exist.")
    } else if(!file.exists(samtools_dir) | !file.exists(bcftools_dir)) {
      print("Indicated samtools and bcftools do not exist. ")
    } else if(file.exists(rnaexp_file) & file.exists(rnabam_file)
       & file.exists(samtools_dir) & file.exists(bcftools_dir)) {
      GenerateListForGetRNASeq(output_peptide_txt_file, width = width)
      output_file_rna_list<-paste(output_peptide_txt_file, ".list.txt", sep="")

      print(paste(samtools_dir, "mpileup -l", output_file_rna_list, "-uf", refdna_file, rnabam_file,
                   ">", paste(output_peptide_txt_file, "list.mp", sep=".")))
      error<-tryCatch2(system(paste(samtools_dir, "mpileup -l", output_file_rna_list, "-uf", refdna_file, rnabam_file,
                   ">", paste(output_peptide_txt_file, "list.mp", sep="."))))

      if(error != 0) skip = TRUE
      if(!skip) {
        print(paste(bcftools_dir, "view -c", paste(output_peptide_txt_file, "list.mp", sep="."),
                    ">", paste(output_peptide_txt_file, "list.vcf", sep=".")))
        error<-tryCatch2(system(paste(bcftools_dir, "view -c", paste(output_peptide_txt_file, "list.mp", sep="."),
                   ">", paste(output_peptide_txt_file, "list.vcf", sep="."))))
      }
    }
  }
  if(indel){
    GetRNAseq_indel(output_peptide_txt_file = output_peptide_txt_file,
                    rnaexp_file = rnaexp_file,
                    output_file_rna_vcf = paste(output_peptide_txt_file, "list.vcf", sep="."))
  } else {
    GetRNAseq(output_peptide_txt_file = output_peptide_txt_file,
              rnaexp_file = rnaexp_file,
              output_file_rna_vcf = paste(output_peptide_txt_file, "list.vcf", sep="."))
  }
}

GetRNAseq_indel<-function(output_peptide_txt_file,
                          rnaexp_file,
                          output_file_rna_vcf){
  data<-t(sapply(scan(output_peptide_txt_file, "character", sep="\n"), function(x) strsplit(x, "\t")[[1]]))
  if(ifelse(is.na(rnaexp_file), TRUE, !file.exists(rnaexp_file))){
    ratio_matrix<-matrix(nrow=nrow(data), ncol=3, NA)
    write.table(cbind(data, ratio_matrix), output_peptide_txt_file,
                row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")
    return()
  }

  ##Get RNA Data
  temp<-scan(rnaexp_file, "character", sep="\n")
  if(length(temp) < 2){
    ratio_matrix<-matrix(nrow=nrow(data), ncol=3, NA)
    write.table(cbind(data, ratio_matrix), commandArgs(TRUE)[1],
                row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")
    return()
  }
  rna<-t(sapply(temp, function(x) strsplit(x, "\t")[[1]]))
  rna<-rna[-1,]
  rna_pos<-t(sapply(rna[,2], function(x) strsplit(x, ":|-")[[1]]))
  hit<-match(sapply(data[,2], function(x) strsplit(x,"_")[[1]][2]), rna[,1])

  #data[,chr], data[,m_position]
  for(x in which(is.na(hit))){
    hit_pos<-which(data[x,3]==rna_pos[,1])
    abs_hit_pos<-hit_pos[which(as.numeric(data[x, 12]) > as.numeric(rna_pos[hit_pos, 2]) &
                               as.numeric(data[x, 12]) < as.numeric(rna_pos[hit_pos, 3]))]
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
  print(output_file_rna_vcf)
  if(ifelse(is.na(output_file_rna_vcf), FALSE, file.exists(output_file_rna_vcf))){
    ratio<-t(sapply(scan(output_file_rna_vcf, "character", sep="\n"),
                    function(x) strsplit(x, "\t")[[1]]))
    ratio<-t(sapply(scan(output_file_rna_vcf, "character", sep="\n", skip=which(lapply(ratio, length)>2)[1]),
                    function(x) strsplit(x, "\t")[[1]]))
    if(ncol(ratio)==0){
      ratio_matrix<-matrix(nrow=nrow(data), ncol=2, NA)
    } else {
      size<-sum(as.numeric(sapply(strsplit(ratio[1,3],"M|D|I")[[1]],
                                  function(x) rev(strsplit(x,"N|S|H|P")[[1]])[1])))
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
  } else {
    ratio_matrix<-matrix(nrow=nrow(data), ncol=2, NA)
  }
  write.table(cbind(data, ratio_matrix), output_peptide_txt_file,
              row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")
}

GetRNAseq<-function(output_peptide_txt_file,
                    rnaexp_file,
                    output_file_rna_vcf){
  data<-t(sapply(scan(output_peptide_txt_file, "character", sep="\n"), function(x) strsplit(x, "\t")[[1]]))

  if(ifelse(is.na(rnaexp_file), TRUE, !file.exists(rnaexp_file))){
    ratio_matrix<-matrix(nrow=nrow(data), ncol=3, NA)
    write.table(cbind(data, ratio_matrix), output_peptide_txt_file,
                row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")
    return()
  }

  ##Get RNA Data
  temp<-scan(rnaexp_file, "character", sep="\n")
  if(length(temp) < 2){
    ratio_matrix<-matrix(nrow=nrow(data), ncol=3, NA)
    write.table(cbind(data, ratio_matrix), commandArgs(TRUE)[1],
                row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")
    return()
  }
  rna<-t(sapply(temp, function(x) strsplit(x, "\t")[[1]]))
  rna<-rna[-1,]
  rna_pos<-t(sapply(rna[,2], function(x) strsplit(x, ":|-")[[1]]))
  hit<-match(sapply(data[,2], function(x) strsplit(x,"_")[[1]][2]), rna[,1])

  #data[,chr], data[,m_position]
  for(x in which(is.na(hit))){
    hit_pos<-which(data[x,3]==rna_pos[,1])
    abs_hit_pos<-hit_pos[which(as.numeric(data[x,12]) > as.numeric(rna_pos[hit_pos,2]) &
                               as.numeric(data[x,12]) < as.numeric(rna_pos[hit_pos,3]))]
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
  print(output_file_rna_vcf)
  if(ifelse(is.na(output_file_rna_vcf), FALSE, file.exists(output_file_rna_vcf))){
    ratio<-t(sapply(scan(output_file_rna_vcf, "character", sep="\n"),
                    function(x) strsplit(x, "\t")[[1]]))
    ratio<-t(sapply(scan(output_file_rna_vcf, "character", sep="\n", skip=which(lapply(ratio, length)>2)[1]),
                    function(x) strsplit(x, "\t")[[1]]))
    if(ncol(ratio)==0){
      ratio_matrix<-matrix(nrow=nrow(data), ncol=2, NA)
    } else {

      ratio_matrix<-NULL
      for(i in 1:nrow(data)){
        hit<-which(data[i, 3]==ratio[,1] & data[i, 12]==ratio[,2])
        if(length(hit)==1 & !is.na(data[i,13])){
          l<-strsplit(ratio[hit,8], "=|,|;")[[1]]
          hit2<-grep("DP4",l)
          if(length(hit2)==1) count<-as.numeric(l[(hit2+1):(hit2+4)])
          else count<-c(.1, 0, 0, 0)
        } else {
          count<-c(.1, 0, 0, 0)
        }
        ratio_matrix<-rbind(ratio_matrix,
                            c(paste(c(sum(count[c(3,4)]),sum(count)), collapse="/"),
                              ifelse(is.na(data[i,tail_col]), 0, as.numeric(data[i, tail_col])) * sum(count[c(3,4)])/sum(count)))
      }
    }
  } else {
    ratio_matrix<-matrix(nrow=nrow(data), ncol=2, NA)
  }
  write.table(cbind(data, ratio_matrix), output_peptide_txt_file,
              row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")
}

GenerateListForGetRNASeq<-function(output_peptide_txt_file,
                                   width){
  data<-t(sapply(scan(output_peptide_txt_file, "character", sep="\n"),
                 function(x) strsplit(x, "\t")[[1]]))
  data<-cbind(data[,3], as.numeric(data[,12]) - width, as.numeric(data[,12]) + width)
  write.table(data, paste(output_peptide_txt_file, ".list.txt", sep=""), row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")
}
