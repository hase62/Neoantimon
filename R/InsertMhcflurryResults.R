InsertMhcflurryResults<-function(result,
                                 output_peptide_prefix,
                                 hla_types){
  print("Writing mhcflurry Results...")

  result_mhcflu <- as.matrix(result)
  for(h_ in 1:length(hla_types)){
    if(file.exists(paste(output_peptide_prefix, ".mhcflurry.HLACLASS1.", h_, ".peptide.mhcflurry.csv", sep = ""))){
      if(requireNamespace("data.table", quietly=TRUE)) {
        mhcf_mut <-data.table::fread(paste(output_peptide_prefix, ".mhcflurry.HLACLASS1.", h_, ".peptide.mhcflurry.csv", sep = ""),
                         stringsAsFactors=FALSE, sep=",", data.table = FALSE)
      } else {
        index <- scan(paste(output_peptide_prefix, ".mhcflurry.HLACLASS1.", h_, ".peptide.mhcflurry.csv", sep = ""),
                      "character", sep = ",", nlines = 1)
        mhcf_mut  <- matrix(scan(paste(output_peptide_prefix, ".mhcflurry.HLACLASS1.", h_, ".peptide.mhcflurry.csv", sep = ""),
                                "character", sep = ",", skip = 1), ncol = length(index), byrow = TRUE)
        colnames(mhcf_mut) <- index
      }
      hla_hit <- which(!is.na(match(result_mhcflu[, 1], mhcf_mut$allele[1])))
      peptide_mut_hit <- match(result_mhcflu[hla_hit, match("Evaluated_Mutant_Peptide", colnames(result_mhcflu))], mhcf_mut$peptide)
      result_mhcflu[hla_hit, match("Mut_EL", colnames(result_mhcflu))] <- as.numeric(mhcf_mut$affinity[peptide_mut_hit])
      result_mhcflu[hla_hit, match("Mut_Rank", colnames(result_mhcflu))] <- mhcf_mut$percentile_rank[peptide_mut_hit]
      if(file.exists(paste(output_peptide_prefix, ".mhcflurry.HLACLASS1.", h_, ".wtpeptide.mhcflurry.csv", sep = ""))){
        if(requireNamespace("data.table", quietly=TRUE)) {
          mhcf_wd <-data.table::fread(paste(output_peptide_prefix, ".mhcflurry.HLACLASS1.", h_, ".wtpeptide.mhcflurry.csv", sep = ""),
                          stringsAsFactors=FALSE, sep=",", data.table = FALSE)
        } else {
          index <- scan(paste(output_peptide_prefix, ".mhcflurry.HLACLASS1.", h_, ".wtpeptide.mhcflurry.csv", sep = ""),
                        "character", sep = ",", nlines = 1)
          mhcf_mut  <- matrix(scan(paste(output_peptide_prefix, ".mhcflurry.HLACLASS1.", h_, ".wtpeptide.mhcflurry.csv", sep = ""),
                                   "character", sep = ",", skip = 1), ncol = length(index), byrow = TRUE)
          colnames(mhcf_mut) <- index
        }
        peptide_wd_hit <-  match(result_mhcflu[hla_hit, match("Evaluated_Wt_Peptide", colnames(result_mhcflu))], mhcf_wd$peptide)
        result_mhcflu[hla_hit, match("Wt_EL",  colnames(result_mhcflu))] <- mhcf_wd$affinity[peptide_wd_hit]
        result_mhcflu[hla_hit, match("Wt_Rank",  colnames(result_mhcflu))] <- mhcf_wd$percentile_rank[peptide_wd_hit]
      }
    } else {
      if(file.exists(paste(output_peptide_prefix, ".HLACLASS1.", h_, ".peptide.txt", sep = ""))){
        tmp <- scan(paste(output_peptide_prefix, ".HLACLASS1.", h_, ".peptide.txt", sep = ""), "character", nlines = 300, sep = "\t")
        hla_hit <- which(!is.na(match(result_mhcflu[, 1], grep("HLA-", strsplit(tmp[length(tmp)], " ")[[1]], value = TRUE))))
        result_mhcflu[hla_hit, match("Mut_EL", colnames(result_mhcflu))] <- NA
        result_mhcflu[hla_hit, match("Mut_Rank", colnames(result_mhcflu))] <- NA
        if(file.exists(paste(output_peptide_prefix, ".HLACLASS1.", h_, ".wtpeptide.txt", sep = ""))){
          result_mhcflu[hla_hit, match("Wt_EL", colnames(result_mhcflu))] <- NA
          result_mhcflu[hla_hit, match("Wt_Rank", colnames(result_mhcflu))] <- NA
        }
      }
    }
  }
  write.table(result_mhcflu, paste(output_peptide_prefix, ".HLACLASS1.ALL.mhcflurry.txt", sep = ""),
              col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
  return(result_mhcflu)
}
