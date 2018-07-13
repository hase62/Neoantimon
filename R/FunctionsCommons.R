CheckRequiredFiles<-function(input_file,
                             hla_file,
                             refflat_file,
                             refmrna_file){
  if(!file.exists(input_file)) {
    print(paste("Did not find Mutation File:", input_file))
    return(TRUE)
  }
  if(!file.exists(hla_file)) {
    print(paste("Did not find HLA Table:", hla_file))
    return(TRUE)
  }
  if(!file.exists(refflat_file)) {
    print(paste("Did not find refFlat File:", refflat_file))
    return(TRUE)
  }
  if(!file.exists(refmrna_file)) {
    print(paste("Did not find refMrna File:", refmrna_file))
    return(TRUE)
  }
  return(FALSE)
}

CheckRequiredFiles2 <- function(input_sequence,
                                input_nm_id,
                                hla_file,
                                refflat_file,
                                refmrna_file,
                                nm_id = NA,
                                gene_symbol = NA,
                                reading_frame = 1){
  if(is.na(input_sequence) & is.na(input_nm_id)) {
    print("Sequence is NaN")
    return(TRUE)
  } else if(!is.na(input_sequence) & !is.na(input_nm_id)) {
    print(paste("Please Specify Either One of", input_sequence, "or", input_nm_id))
    return(TRUE)
  } else if(is.na(input_sequence)) {
    print(paste("Sequence:", input_sequence))
  } else if(is.na(input_nm_id)) {
    print(paste("NM_ID:", input_nm_id))
  }
  if(!file.exists(hla_file)) {
    print(paste("Did not find HLA Table:", hla_file))
    return(TRUE)
  } else {
    print(paste("HLA file:", hla_file))
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
  if(length(nm_id) <= 1 & length(gene_symbol) <= 1){
    if(is.na(nm_id) & is.na(gene_symbol)) {
      print("Please Specify Either One of: nm_id or gene_symbol.")
      return(TRUE)
    }else {
      print(paste("NM_ID: ", nm_id))
      print(paste("Gene Symbol:", gene_symbol))
    }
  } else {
    print(paste("NM_ID: ", nm_id))
    print(paste("Gene Symbol:", gene_symbol))
  }
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
  index<-scan(input_file, "character", nlines = 1)

  if(is.na(chr_column)) {
    chr_column<-grep("chr", tolower(index))[1];
    if(is.na(chr_column)) {
      print("Please Manually Indicate chr_column")
      return(0)
    }
  }
  print(paste("Set chr_column as", chr_column))

  if(is.na(mutation_start_column)) {
    mutation_start_column<-grep("start", tolower(index))[1];
    if(is.na(mutation_start_column)) {
      print("Please Manually Indicate mutation_start_column")
      return(0)
    }
  }
  print(paste("Set mutation_start_column as", mutation_start_column))

  if(is.na(mutation_end_column)) {
    mutation_end_column<-grep("end", tolower(index))[1];
    if(is.na(mutation_end_column)) {
      print("Please Manually Indicate mutation_end_column")
      return(0)
    }
  }
  print(paste("Set mutation_end_column as", mutation_end_column))

  if(is.na(mutation_ref_column)) {
    mutation_ref_column<-grep("ref", tolower(index))[1];
    if(is.na(mutation_ref_column)) {
      print("Please Manually Indicate mutation_ref_column")
      return(0)
    }
  }
  print(paste("Set mutation_ref_column as", mutation_ref_column))

  if(is.na(mutation_alt_column)) {
    mutation_alt_column<-grep("alt", tolower(index))[1];
    if(is.na(mutation_alt_column)) {
      print("Please Manually Indicate mutation_alt_column")
      return(0)
    }
  }
  print(paste("Set mutation_alt_column as", mutation_alt_column))

  if(is.na(depth_normal_column)) {
    depth_normal_column<-intersect(grep("depth", tolower(index)), grep("normal", tolower(index)))[1];
    if(is.na(depth_normal_column)) {
      print("Please Manually Indicate depth_normal_column")
    }
  }
  print(paste("Set depth_normal_column as", depth_normal_column))

  if(is.na(depth_tumor_column)) {
    depth_tumor_column<-intersect(grep("depth", tolower(index)), grep("tumor", tolower(index)))[1];
    if(is.na(depth_tumor_column)) {
      print("Please Manually Indicate depth_tumor_column")
    }
  }
  print(paste("Set depth_tumor_column as", depth_tumor_column))

  if(is.na(nm_id_column)) {
    nm_id_column<-grep("AAChange.refGene", tolower(index))[1];
    if(is.na(nm_id_column)) {
      index<-scan(input_file, "character", nlines = 1, skip = 1)
      nm_id_column<-grep("nm_|nr_", tolower(index))[1];
      if(is.na(nm_id_column)) {
        print("Please Manually Indicate nm_id_column")
        return(0)
      }
    }
  }
  print(paste("Set nm_id_column as", nm_id_column))

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


