Read_files <- function(input_annovar_format_file, input_vep_format_file, input_vcf_format_file_and_vep){
  num_files <- length(which(is.na(c(input_annovar_format_file,
                                     input_vep_format_file,
                                     input_vcf_format_file_and_vep))))
  if(num_files != 2) {
    print("Please indicate one of input_annovar_format_file, input_vep_format_file, and input_vcf_format_file_and_vep. ")
    return(TRUE)
  }
  return(FALSE)
}

CheckRequiredFiles<-function(input_file,
                             hla_types,
                             refflat_file,
                             refmrna_file){
  if(is.list(input_file) | is.matrix(input_file)){
    print("Input Data is")
    print(head(input_file))
  } else if(!file.exists(input_file)) {
    print(paste("Did not find Mutation File:", input_file))
    return(TRUE)
  }

  if(is.na(hla_types[1])) {
    print(paste("Did not find HLA Types"))
    return(TRUE)
  } else {
    print(paste("HLAtype:", hla_types))
  }

  if(is.list(refflat_file) | is.matrix(refflat_file)){
    print("refFlat is")
    print(head(refflat_file))
  } else if(!file.exists(refflat_file)) {
    print(paste("Did not find refFlat File:", refflat_file))
    return(TRUE)
  }

  if(is.list(refmrna_file) | is.matrix(refmrna_file)){
    print("refMrna is")
    print(head(input_file))
  } else if(!file.exists(refmrna_file)) {
    print(paste("Did not find refMrna File:", refmrna_file))
    return(TRUE)
  }
  return(FALSE)
}

CheckRequiredFiles2 <- function(input_sequence,
                                input_nm_id,
                                hla_types,
                                refflat_file,
                                refmrna_file,
                                reference_nm_id = NA,
                                reference_gene_symbol = NA,
                                reading_frame = 1){
  if(is.na(input_sequence[1]) & is.na(input_nm_id[1])) {
    print("Sequence is NaN")
    return(TRUE)
  }
  if(is.na(input_sequence[1])) {
    print(paste("NM_ID:", input_nm_id[1]))
  } else if(is.na(input_nm_id[1])) {
    print(paste("Sequence:", input_sequence))
  }
  if(is.na(hla_types[1])) {
    print(paste("Did not find HLA Types"))
    return(TRUE)
  } else {
    print(paste("HLAtype:", hla_types))
  }

  if(!file.exists(refflat_file)) {
    print(paste("Did not find refFlat File:", refflat_file))
    return(TRUE)
  } else {
    print(paste("refFLAT:", refflat_file))
  }
  if(!file.exists(refmrna_file)) {
    print(paste("Did not find refMrna File:", refmrna_file))
    return(TRUE)
  } else {
    print(paste("refMrna:", refmrna_file))
  }
  print(paste("Wt-NM_ID: ", reference_nm_id))
  print(paste("Wt-Gene Symbol:", reference_gene_symbol))

  print(paste("Start Position of Reading Frame is:", reading_frame))
  return(FALSE)
}

CheckRequiredColumns<-function(input_file,
                               chr_column,
                               mutation_start_column,
                               mutation_end_column,
                               mutation_ref_column,
                               mutation_alt_column,
                               nm_id_column,
                               depth_normal_column,
                               depth_tumor_column
                               ){
  if(is.list(input_file) | is.matrix(input_file)){
    index <- colnames(input_file)
  } else {
    index<-scan(input_file, "character", nlines = 1)
  }

  if(is.na(chr_column)) {
    chr_column<-grep("chr", index, ignore.case = TRUE)[1];
    if(is.na(chr_column)) {
      print("Please Manually Indicate chr_column")
      return(0)
    }
  }
  print(paste("Set chr_column as", chr_column, index[chr_column]))

  if(is.na(mutation_start_column)) {
    mutation_start_column<-grep("start", index, ignore.case = TRUE)[1];
    if(is.na(mutation_start_column)) {
      print("Please Manually Indicate mutation_start_column")
      return(0)
    }
  }
  print(paste("Set mutation_start_column as", mutation_start_column, index[mutation_start_column]))

  if(is.na(mutation_end_column)) {
    mutation_end_column<-grep("end", index, ignore.case = TRUE)[1];
    if(is.na(mutation_end_column)) {
      print("Please Manually Indicate mutation_end_column")
      return(0)
    }
  }
  print(paste("Set mutation_end_column as", mutation_end_column, index[mutation_end_column]))

  if(is.na(mutation_ref_column)) {
    mutation_ref_column<-grep("ref", index, ignore.case = TRUE)[1];
    if(is.na(mutation_ref_column)) {
      print("Please Manually Indicate mutation_ref_column")
      return(0)
    }
  }
  print(paste("Set mutation_ref_column as", mutation_ref_column, index[mutation_ref_column]))

  if(is.na(mutation_alt_column)) {
    mutation_alt_column<-grep("alt", index, ignore.case = TRUE)[1];
    if(is.na(mutation_alt_column)) {
      print("Please Manually Indicate mutation_alt_column")
      return(0)
    }
  }
  print(paste("Set mutation_alt_column as", mutation_alt_column, index[mutation_alt_column]))

  if(is.na(depth_normal_column)) {
    depth_normal_column<-intersect(grep("depth", index, ignore.case = TRUE),
                                   grep("normal", index, ignore.case = TRUE))[1];
    if(is.na(depth_normal_column)) {
      print("Please Manually Indicate depth_normal_column if you want to use.")
    }
  }
  print(paste("Set depth_normal_column as", depth_normal_column, index[depth_normal_column]))

  if(is.na(depth_tumor_column)) {
    depth_tumor_column<-intersect(grep("depth", index, ignore.case = TRUE),
                                  grep("tumor", index, ignore.case = TRUE))[1];
    if(is.na(depth_tumor_column)) {
      print("Please Manually Indicate depth_tumor_column if you want to use.")
    }
  }
  print(paste("Set depth_tumor_column as", depth_tumor_column, index[depth_tumor_column]))

  if(is.na(nm_id_column)) {
    nm_id_column<-grep("AAChange", index, ignore.case = TRUE)[1];
    if(is.na(nm_id_column)) {
      index<-scan(input_file, "character", nlines = 1, skip = 1)
      nm_id_column<-grep("nm_|nr_", tolower(index))[1];
      if(is.na(nm_id_column)) {
        print("Please Manually Indicate nm_id_column")
        return(0)
      }
    }
  }
  print(paste("Set nm_id_column as", nm_id_column, index[nm_id_column]))

  return(c(chr_column,
           mutation_start_column,
           mutation_end_column,
           mutation_ref_column,
           mutation_alt_column,
           nm_id_column,
           depth_normal_column,
           depth_tumor_column))
}

SettingNetMHCpan<-function(netMHCpan_dir){
  netMHCpan_script<-scan(netMHCpan_dir, "character", sep="\n", blank.lines.skip = FALSE)
  netMHCpan_par<-gsub("\\./", "/", paste(rev(rev(strsplit(netMHCpan_dir, "/")[[1]])[-1]), collapse = "/"))
  if(!file.exists(paste(getwd(), "/", netMHCpan_par, "/tmp", sep=""))){
    dir.create(paste(getwd(), "/", netMHCpan_par, "/tmp", sep=""))
  }
  netMHCpan_script<-gsub("#setnv", "setenv", netMHCpan_script)

  netMHCpan_script[grep("setenv\tNMHOME", netMHCpan_script)]<-
    paste("setenv\tNMHOME ", getwd(), "/", netMHCpan_par, sep="")
  write.table(netMHCpan_script, netMHCpan_dir, row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")
}

SettingNetMHCIIpan<-function(netMHCIIpan_dir){
  netMHCpan_script<-scan(netMHCIIpan_dir, "character", sep="\n", blank.lines.skip = FALSE)
  netMHCpan_par<-gsub("\\./","/", paste(rev(rev(strsplit(netMHCIIpan_dir, "/")[[1]])[-1]), collapse = "/"))
  if(!file.exists(paste(getwd(), "/", netMHCpan_par, "/tmp", sep=""))){
    dir.create(paste(getwd(), "/", netMHCpan_par, "/tmp", sep=""))
  }
  netMHCpan_script[grep("setenv\tNMHOME", netMHCpan_script)]<-
    paste("setenv\tNMHOME ", getwd(), "/", netMHCpan_par, sep="")
  netMHCpan_script[grep("setenv\tNMHOME", netMHCpan_script) + 1]<-
    paste("setenv\tTMPDIR ", getwd(), "/", netMHCpan_par, "/tmp", sep="")
  write.table(netMHCpan_script, netMHCIIpan_dir, row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")
}

tryCatch2<-function(g){
  tryCatch(g,
           error=function(e){
             message(e)
             cat("sasasa")
           }
  )
}

apply2<-function(x, y, f){
  if(is.null(nrow(x))){
    return(f(x))
  } else {
    apply(x, y, function(a) f(a))
  }
}

apply3<-function(x, y, z){
  if(is.null(dim(x))){
    return(x)
  } else {
    return(apply(x, y, z))
  }
}

apply4<-function(x, y, z){
  out <- apply(x, y, z)
  if(!is.matrix(out)) out <- t(out)
  return(out)
}

cut_peptide <- function(peptide, lens) {
  pep_result <- NULL
  for(l in lens){
    tmp <- sapply(1:(nchar(peptide) - l + 1), function(x) substr(peptide, x, x + l - 1))
    pep_result <- rbind(pep_result, cbind(1:length(tmp), tmp))
  }
  return(pep_result)
}

getHLAtypes<-function(hla_file, file_name_in_hla_table){
  hla <- t(sapply(scan(hla_file, "character", sep="\n"), function(x) strsplit(x, "\t")[[1]]))
  hit <- match(file_name_in_hla_table, hla[, 1])
  if(is.na(hit)) {
    print(file_name_in_hla_table, "is not included in", hla_file)
    return (NA)
  }
  return(hla[hit, -1])
}

check_dna_validity <- function(dna, nm_id, exon_end, exon_start, ambiguous_between_exon, Last, Pass){
  if(is.na(dna)){
    print(paste(nm_id, "was not found in refMrn."))
    return(TRUE)
  }
  if(nchar(dna) != sum(exon_end - exon_start)){
    dif<-sum(exon_end - exon_start) - nchar(dna)
    if(abs(dif) <= ambiguous_between_exon) {
      if(Last && !Pass){
        print(paste("cDNA Length does not Match to Exon-Start/End Length, Skip", nm_id))
      } else {
        print(paste("Permit Ambiguous Exonic Region:", nm_id))
      }
      return(TRUE)
    }
  }
  return(FALSE)
}

get_relative_mutation_position <- function(strand, exon_end, m_start, exon_start){
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
  return(m_point)
}

get_relative_translation_start_position <- function(strand, exon_end, trans_start, exon_start, trans_end){
  if(strand=="+"){
    point<-(exon_end > trans_start)
    ts_point<-sum((exon_end - exon_start)[!point]) +
      (trans_start - exon_start[which(point)[1]]) + 1
  }else{
    point<-(exon_start > trans_end)
    ts_point<-sum((exon_end - exon_start)[point]) +
      (exon_end[rev(which(!point))[1]] - trans_end) + 1
  }
  return(ts_point)
}

check_start_codon <- function(dna, ts_point, ambiguous_codon, nm_id){
  d<-0
  if(substr(dna, ts_point, ts_point + 2) != "atg"){
    flag <- FALSE
    for(d in  (-1 * ambiguous_codon):(ambiguous_codon)){
      if(substr(dna, ts_point + d, ts_point + 2 + d)=="atg"){
        flag<-TRUE
        return(NULL)
      }
    }
    if(flag){
      if(d < 0){
        dna <- sub(" ", "", paste(paste(rep("x", -d), collapse = ""), dna, collapse = ""))
      }else{
        dna <- substr(dna, d + 1, nchar(dna))
      }
      print("Permit Ambiguous Codon Start")
    }else{
      print(paste("Start Position is not ATG, Skip", nm_id))
      return(-999)
    }
  }
  return(d)
}

get_relative_translation_end_position <- function(strand, exon_end, trans_start, exon_start, trans_end){
  if(strand == "+"){
    point<-(exon_end >= trans_end)
    te_point<-sum((exon_end - exon_start)[!point]) +
      (trans_end - exon_start[which(point)[1]])
  }else{
    point<-(exon_start > trans_start)
    te_point<-sum((exon_end - exon_start)[point]) +
      (exon_end[rev(which(!point))[1]] - trans_start)
  }
  return(te_point)
}

check_stop_codon <- function(dna, te_point, ts_point, ambiguous_codon, amino, nm_id){
  e<-0
  if(amino[match(substr(dna, te_point-2, te_point), codon)]!="X"){
    dna_trans<-substr(dna, ts_point, te_point)
    flag <- FALSE
    dif <- nchar(dna_trans)%%3
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
      return(-999)
    } else {
      print("Permit Ambiguous Codon Start")
    }
  }
  return(e)
}

make_normal_peptide <- function(dna_trans, amino, codon, k, e){
  peptide_normal <- NULL
  while(nchar(dna_trans) >= 3){
    peptide_normal <- c(peptide_normal, amino[match(substr(dna_trans, 1, 3), codon)])
    dna_trans <- substr(dna_trans, 4, nchar(dna_trans))
  }
  if(k==e & match("X", peptide_normal) < length(peptide_normal)){
    return(NULL)
  }
  return(peptide_normal)
}

make_mutated_dna <- function(strand, dna_trans, m_point_2, m_ref, m_alt, trans_to, trans_from){
  if(strand=="+"){
    if(substr(dna_trans, m_point_2, m_point_2) == tolower(m_ref))
      substr(dna_trans, m_point_2, m_point_2) <- tolower(m_alt)
  } else {
    if(substr(dna_trans, m_point_2, m_point_2) == trans_to[match(tolower(m_ref), trans_from)])
      substr(dna_trans, m_point_2, m_point_2) <- trans_to[match(tolower(m_alt), trans_from)]
  }
  return(dna_trans)
}

make_mutated_peptide <- function(dna_trans_mut, amino, codon){
  peptide <- NULL
  while(nchar(dna_trans_mut) >= 3){
    a <- amino[match(substr(dna_trans_mut, 1, 3), codon)]
    peptide <- c(peptide, a)
    if(a=="X") break
    dna_trans_mut <- substr(dna_trans_mut, 4, nchar(dna_trans_mut))
  }
  return(peptide)
}

generate_fraction <- function(m_point_2, max_peptide_length, peptide){
  peptide_start <- ceiling(m_point_2 / 3.0) - max_peptide_length
  if(peptide_start < 1) peptide_start <- 1
  peptide_end <- ceiling(m_point_2 / 3.0) + max_peptide_length
  if(peptide_end > length(peptide)) peptide_end <- length(peptide)
  return(peptide_start:peptide_end)
}

generate_fraction_indel <- function (peptide, peptide_normal, max_peptide_length){
  min_len <- min(length(peptide), length(peptide_normal))
  peptide_start <- which(peptide[1:min_len] != peptide_normal[1:min_len])[1] - max_peptide_length + 1
  if(is.na(peptide_start)) {
    print("Peptide Start is NA")
    return(NULL)
  }
  if(peptide_start < 1) peptide_start<-1
  peptide_end<-which(rev(peptide)[1:min_len] != rev(peptide_normal)[1:min_len])[1]
  if(is.na(peptide_end)) peptide_end <- min_len
  peptide_end <- min_len - peptide_end + max_peptide_length + 10
  if(peptide_end > length(peptide)) peptide_end = length(peptide)
  peptide <- peptide[peptide_start:peptide_end]
  peptide_end <- peptide_end + 10
  if(peptide_end > length(peptide_normal)) peptide_end = length(peptide_normal)
  peptide_normal <- peptide_normal[peptide_start:min(peptide_end, length(peptide_normal))]
  return(list(peptide, peptide_normal))
}

integrate_same_peptide <- function(refFasta, fasta, fasta_wt){
  i <- 1
  while(i <= nrow(refFasta)){
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
    temp1<-paste(refFasta[hit,3], collapse=";")
    temp2<-paste(refFasta[hit,4], collapse=";")
    temp3<-paste(refFasta[hit,9], collapse=";")
    temp4<-paste(refFasta[hit,10], collapse=";")
    refFasta[i,3]<-temp1
    refFasta[i,4]<-temp2
    refFasta[i,9]<-temp3
    refFasta[i,10]<-temp4
    refFasta<-refFasta[-hit[-1],]
    if(is.null(nrow(refFasta))){
      refFasta<-t(refFasta)
    }
    fasta<-fasta[-(c(hit[-1] * 2 - 1, hit[-1] * 2))]
    fasta_wt<-fasta_wt[-(c(hit[-1] * 2 - 1, hit[-1] * 2))]
    i<-i+1
  }
  return(list(refFasta, fasta, fasta_wt))
}

check_multiple_snvs <- function(data, multiple_variants, i, exon_start, mutation_start_column, chr_column, exon_end, chr){
  multi_i <- integer(0)
  if(multiple_variants & nrow(data) > 1){
    multi_i <- (1:nrow(data))[-i][sapply((1:nrow(data))[-i],
                                         function(x) length(which(exon_start < data[x, mutation_start_column] &
                                                                    data[x, mutation_start_column] <= exon_end  &
                                                                    data[x, chr_column] == chr)) == 1)]
  }
  return(multi_i)
}

apply_multiple_snvs <- function(data, multiple_variants, i, exon_start, mutation_start_column, chr_column, mutation_ref_column, mutation_alt_column, exon_end, chr, strand, dna_trans_mut, trans_to, trans_from){
  multi_i <- check_multiple_snvs(data, multiple_variants, i, exon_start, mutation_start_column, chr_column, exon_end, chr)
  for(multi_i_element in multi_i){
    m_point_2_mv <- get_relative_mutation_position(strand, exon_end, as.numeric(data[multi_i_element, mutation_start_column]), exon_start)
    dna_trans_mut <- make_mutated_dna(strand, dna_trans_mut, m_point_2_mv, data[multi_i_element, mutation_ref_column], data[multi_i_element, mutation_alt_column], trans_to, trans_from)
  }
  return(dna_trans_mut)
}

check_multiple_snps <- function(SNPs_vcf, exon_start, mutation_start_column, exon_end, chr){
  multi_i <- integer(0)
  if(nrow(SNPs_vcf) > 1){
    multi_i <- (1:nrow(SNPs_vcf))[sapply((1:nrow(SNPs_vcf)),
                                         function(x) length(which(exon_start < SNPs_vcf[x, 2] &
                                                                  SNPs_vcf[x, 2] <= exon_end  &
                                                                  SNPs_vcf[x, 1] == chr)) == 1)]
  }
  return(multi_i)
}

apply_multiple_snps <- function(SNPs_vcf, exon_start, mutation_start_column, exon_end, chr, strand, dna_trans, trans_to, trans_from){
  multi_i <- check_multiple_snps(SNPs_vcf, exon_start, mutation_start_column, exon_end, chr)
  for(multi_i_element in multi_i){
    m_point_2_mv <- get_relative_mutation_position(strand, exon_end, as.numeric(SNPs_vcf[multi_i_element, 2]), exon_start)
    dna_trans <- make_mutated_dna(strand, dna_trans, m_point_2_mv, tolower(SNPs_vcf[multi_i_element, 4]), tolower(SNPs_vcf[multi_i_element, 5]), trans_to, trans_from)
  }
  return(dna_trans)
}

make_indel_dna <- function(strand, dna_trans, m_point_2, m_alt, trans_to, trans_from, m_ref){
  if(strand == "+"){
    if(m_ref == "-"){
      #Insertion
      dna_trans<-paste(substr(dna_trans, 1, m_point_2 - 1),
                       paste(sapply(substring(m_alt, 1:nchar(m_alt), 1:nchar(m_alt)),
                                    function(x) trans_to[match(tolower(x), trans_from)]), collapse=""),
                       substr(dna_trans, m_point_2, nchar(dna_trans)), sep="")
    } else {
      ref_ <- substr(dna_trans, m_point_2, m_point_2 + nchar(m_ref) - 1)
      alt_ <- paste(substring(tolower(m_ref), 1:nchar(m_ref), 1:nchar(m_ref)), collapse="")
      match_ <- ref_ == alt_
      match_length <- length(which(strsplit(ref_, "")[[1]] == strsplit(alt_, "")[[1]]))
      match_pval <- 1 - pbinom(match_length, nchar(ref_), 0.25)
      if(!match_ & match_pval < 0.03){
        print(paste("Ref is", ref_, ", and vcf is", alt_))
        print(paste("The Ref and vcf have miss-match ... but p-val(", match_pval, ") is less tnan 0.03 ... Continue Convertion. ", sep = ""))
        match_ <- TRUE
      }
      if(match_){
        dna_trans<-paste(substr(dna_trans, 1, m_point_2 - 1),
                         gsub("NA", "", paste(sapply(substring(m_alt, 1:nchar(m_alt), 1:nchar(m_alt)),
                                                     function(x) trans_to[match(tolower(x),trans_from)]), collapse="")),
                         substr(dna_trans, m_point_2 + nchar(m_ref), nchar(dna_trans)), sep="")
      } else {
        print("The Ref and vcf are not matched")
        return(NULL)
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
      ref_ <- paste(substring(substr(dna_trans, m_point_2 - nchar(m_ref) + 1, m_point_2), 1:nchar(m_ref), 1:nchar(m_ref)), collapse="")
      alt_ <- paste(sapply(rev(substring(tolower(m_ref), 1:nchar(m_ref), 1:nchar(m_ref))), function(x) trans_to[match(tolower(x),trans_from)]), collapse="")
      match_ <- ref_ == alt_
      match_length <- length(which(strsplit(ref_, "")[[1]] == strsplit(alt_, "")[[1]]))
      match_pval <- 1 - pbinom(match_length, nchar(ref_), 0.25)
      if(!match_ & match_pval < 0.03){
        print(paste("Ref is", ref_, ", and vcf is", alt_))
        print(paste("The Ref and vcf have miss-match ... but p-val(", match_pval, ") is less tnan 0.03 ... Continue Convertion. ", sep = ""))
        match_ <- TRUE
      }
      if(match_){
        dna_trans<-paste(substr(dna_trans, 1, m_point_2 - nchar(m_ref)),
                         gsub("NA","",paste(sapply(rev(substring(m_alt, 1:nchar(m_alt), 1:nchar(m_alt))),
                                                   function(x) trans_to[match(tolower(x),trans_from)]), collapse="")),
                         substr(dna_trans, m_point_2 + 1, nchar(dna_trans)), sep="")
      } else {
        print("The Ref and vcf are not matched")
        return(NULL)
      }
    }
  }
  return(dna_trans)
}

check_multiple_snvs_to_indel <- function(data, multiple_variants, exon_start, mutation_start_column, chr_column, exon_end, chr){
  multi_i <- integer(0)
  if(multiple_variants & nrow(data) > 1){
    multi_i <- (1:nrow(data))[sapply((1:nrow(data)),
                                     function(x) length(which(exon_start < data[x, mutation_start_column] &
                                                              data[x, mutation_start_column] <= exon_end  &
                                                              data[x, chr_column] == chr)) == 1)]
  }
  return(multi_i)
}

apply_multiple_snvs_to_indel <- function(data, multiple_variants, exon_start, mutation_start_column, chr_column, mutation_ref_column, mutation_alt_column, exon_end, chr, strand, dna_trans_mut, trans_to, trans_from){
  multi_i <- check_multiple_snvs_to_indel(data, multiple_variants, exon_start, mutation_start_column, chr_column, exon_end, chr)
  for(multi_i_element in multi_i){
    m_point_2_mv <- get_relative_mutation_position(strand, exon_end, as.numeric(data[multi_i_element, mutation_start_column]), exon_start)
    dna_trans_mut <- make_mutated_dna(strand, dna_trans_mut, m_point_2_mv, data[multi_i_element, mutation_ref_column], data[multi_i_element, mutation_alt_column], trans_to, trans_from)
  }
  return(dna_trans_mut)
}

check_valid_indel <- function(peptide, ignore_short, max_peptide_length){
  X <- grep("X", peptide)
  if(length(X) > 0 & ignore_short){
    if(X < 8) {
      print("Indel is Too Short")
      return(TRUE)
    }
  }
  if(max_peptide_length >= 15 & length(X) > 0 & ignore_short){
    if(X < 15) {
      print("Indel is Too Short")
      return(TRUE)
    }
  }
  return(FALSE)
}

read_data <- function(input_file){
  data <- NULL
  tmp <- scan(input_file, "character", sep = "\n", nlines = 200)
  read_start <- grep("\\#chr", tolower(tmp))[1]
  if(is.na(read_start)) read_start <- grep("\\#uploaded_variation", tolower(tmp))[1]
  if(is.na(read_start)) read_start <- rev(grep("\\#", tmp))[1] + 1
  if(is.na(read_start)) read_start <- grep("chr", tolower(tmp))[1]
  if(is.na(read_start)) read_start <- grep("uploaded_variation", tolower(tmp))[1]
  if(is.na(read_start)) read_start <- 1
  print(paste("Please Confirm that Reading Start Line is", read_start))
  print(tmp[read_start])
  if(requireNamespace("data.table", quietly=TRUE)) {
    data <- data.table::fread(input_file, stringsAsFactors=FALSE, sep="\t", skip = read_start - 1, header =TRUE, data.table = FALSE)
  } else {
    index <- scan(input_file, "character", sep = "\t", nlines = 1, skip = read_start - 1)
    data  <- matrix(scan(input_file, "character", sep = "\t", skip = read_start), ncol = length(index), byrow = TRUE)
    colnames(data) <- index
  }
  return(data)
}

read_refFlat <- function(refflat_file){
  if(requireNamespace("data.table", quietly=TRUE)) {
    list_nm <- data.table::fread(refflat_file, stringsAsFactors=FALSE, header = FALSE, sep="\t", data.table = FALSE)
  } else {
    index <- scan(refflat_file, "character", sep = "\t", nlines = 1)
    list_nm  <- matrix(scan(refflat_file, "character", sep = "\t", skip = 1), ncol = length(index), byrow = TRUE)
  }
  return(list_nm)
}

read_refmrn <- function(refmrna_file){
  if(requireNamespace("data.table", quietly=TRUE)) {
    list_mra <- data.table::fread(refmrna_file, stringsAsFactors=FALSE, header = FALSE, sep='\t', data.table = FALSE)[, 1]
  } else {
    list_mra <- scan(refmrna_file, "character", sep = "\n")
  }
  return(list_mra)
}

read_1col_by_fread_or_scan <- function(f_name) {
  if(requireNamespace("data.table", quietly=TRUE)) {
    tmp <- data.table::fread(f_name, stringsAsFactors=FALSE, sep='\n', data.table = FALSE)[, 1]
  } else {
    tmp <- scan(f_name, "character", sep = "\n")
  }
  return(tmp)

}

#library("ensemblVEP")
#param <- VEPFlags()
#f_path <- "data/sample_vcf.snps.vcf"
#file <- system.file(f_path, "ex2.vcf", package="VariantAnnotation")
#gr <- ensemblVEP(file)

