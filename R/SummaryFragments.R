#'Export Summary Count from Indel/SV Results
#'
#'@param Input Input file generated from MainSNVClass1,2.
#'
#'@param mut_IC50_th The threshold for mutant peptide to be neoantigen.
#'
#'@param Total_RNA_th The total RNA expression threshold.
#'
#'@param Tumor_RNA_th The tumor specific RNA expression threshold.
#'
#'@param MutRatio_th The mutation ratio threshold.
#'
#'@param WriteLongIndel If setting a file name, Write Long Indels of which the p-value is less than 0.05.
#'
#'@return Num_Alteration The number of evaluated alterations.
#'
#'@return Num_Alteration_Generating_NeoAg The number of evaluated alterations that can generate neoantigen.
#'
#'@return Num_Peptide The number of evaluated peptifdes.
#'
#'@return Num_Peptide_Generating_NeoAg The number of evaluated peptides that can be neoantigen.
#'
#'@export
Export_Summary_Fragments <- function(Input,
                                     mut_IC50_th,
                                     Total_RNA_th = NA,
                                     Tumor_RNA_th = NA,
                                     MutRatio_th = NA,
                                     WriteLongIndel = NA){

  #Calculate Pvalues
  theoretical_pro <- 3 * 4^3 / 61^2
  pvalues <- sapply(Input[, match("Mutant_Peptide", colnames(Input))],
                    function(x) (1 - theoretical_pro)^nchar(x))
  Input <- cbind(Input, pvalues)
  colnames(Input)[ncol(Input)] <- "Pvalue"

  index <- colnames(Input)

  # Write Long Indels
  if(!is.na(WriteLongIndel)){
    unq <- unique(Input[, match("Mutant_Peptide", index)])
    hit <- unq[as.numeric(Input[match(unq, Input[, match("Mutant_Peptide", index)]), match("Pvalue", index)]) < 0.001]
    if(length(hit) > 0){
      tmp <- Input[match(hit, Input[, match("Mutant_Peptide", index)]), (match("Chr", colnames(Input))):ncol(Input)]
      if(length(hit)==1) tmp <- t(tmp)
      write.table(tmp,
                  WriteLongIndel,
                  row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
    }
  }

  # Count All
  unq_grp <- unique(Input[, match("GroupID", index)])
  unq_nm <- unique(Input[, match("NM_ID", index)])
  unq_pep <- unique(Input[, match("Mutant_Peptide", index)])
  names(unq_grp) <- sapply(unq_grp, function(x) min(as.numeric(Input[which(!is.na(match(Input[, match("GroupID", index)], x))), match("Pvalue", index)])))
  names(unq_nm) <- sapply(unq_nm, function(x) min(as.numeric(Input[which(!is.na(match(Input[, match("NM_ID", index)], x))), match("Pvalue", index)])))
  names(unq_pep) <- sapply(unq_pep, function(x) min(as.numeric(Input[which(!is.na(match(Input[, match("Mutant_Peptide", index)], x))), match("Pvalue", index)])))

  Num_Grp <- length(unq_grp)
  Num_NM <- length(unq_nm)
  Num_Peptide <- length(unq_pep)
  Num_Peptide_Per_Grp <- sapply(unq_grp, function(x) length(which(!is.na(match(Input[, match("GroupID", index)], x)))))
  Num_Peptide_Per_NM <- sapply(unq_nm, function(x) length(which(!is.na(match(Input[, match("NM_ID", index)], x)))))
  Num_Peptide_Per_Pep <- sapply(unq_pep, function(x) length(which(!is.na(match(Input[, match("Mutant_Peptide", index)], x)))))

  # Conditioning
  if(!is.na(Total_RNA_th)){
    Input <- Input[as.numeric(Input[, match("Total_RNA", index)]) > Total_RNA_th, ]
  }
  if(!is.na(Tumor_RNA_th)){
    Input <- Input[as.numeric(Input[, match("Tumor_RNA", index)]) > Tumor_RNA_th, ]
  }
  if(!is.na(MutRatio_th)){
    Input <- Input[as.numeric(Input[, match("MutRatio", index)]) > MutRatio_th, ]
  }

  # Count Conditioned
  Num_Cond_Grp <- length(unq_grp)
  Num_Cond_NM <- length(unq_nm)
  Num_Cond_Peptide <- length(unique(Input[, match("Evaluated_Mutant_Peptide", index)]))
  Num_Cond_Peptide_Per_Grp <- sapply(unq_grp, function(x) length(which(!is.na(match(Input[, match("GroupID", index)], x)))))
  Num_Cond_Peptide_Per_NM <- sapply(unq_nm, function(x) length(which(!is.na(match(Input[, match("NM_ID", index)], x)))))
  Num_Cond_Peptide_Per_Pep <- sapply(unq_pep, function(x) length(which(!is.na(match(Input[, match("Mutant_Peptide", index)], x)))))

  # Extract by IC50
  Input <- Input[as.numeric(Input[, match("Mut_IC50", index)]) < mut_IC50_th, ]

  if(is.null(Input) | is.null(dim(Input)[1])) {
    if(length(Input) > 10) {
      Input <- t(Input)
    }
  }

  # Count Rest
  Num_Rest_Grp <- length(unq_grp)
  Num_Rest_NM <- length(unq_nm)
  Num_Rest_Peptide <- length(unique(Input[, match("Evaluated_Mutant_Peptide", index)]))
  Num_Rest_Peptide_Per_Grp <- sapply(unq_grp, function(x) length(which(!is.na(match(Input[, match("GroupID", index)], x)))))
  Num_Rest_Peptide_Per_NM <- sapply(unq_nm, function(x) length(which(!is.na(match(Input[, match("NM_ID", index)], x)))))
  Num_Rest_Peptide_Per_Pep <- sapply(unq_pep, function(x) length(which(!is.na(match(Input[, match("Mutant_Peptide", index)], x)))))

  ans <- c(Num_Grp, Num_Cond_Grp, Num_Rest_Grp,
           Num_NM, Num_Cond_NM, Num_Rest_NM,
           Num_Peptide, Num_Cond_Peptide, Num_Rest_Peptide)
  names(ans) <- c("Num_Grp", "Num_Cond_Grp", "Num_Rest_Grp",
                  "Num_NM", "Num_Cond_NM", "Num_Rest_NM",
                  "Num_Peptide", "Num_Cond_Peptide", "Num_Rest_Peptide")

  ratio_pep_grp <- round(rbind(Num_Peptide_Per_Grp,
                               Num_Cond_Peptide_Per_Grp,
                               Num_Rest_Peptide_Per_Grp,
                               Num_Rest_Peptide_Per_Grp / Num_Cond_Peptide_Per_Grp,
                               -log10(as.numeric(names(unq_grp)))), 3)

  ratio_pep_nm <- round(rbind(Num_Peptide_Per_NM,
                              Num_Cond_Peptide_Per_NM,
                              Num_Rest_Peptide_Per_NM,
                              Num_Rest_Peptide_Per_NM / Num_Cond_Peptide_Per_NM,
                              -log10(as.numeric(names(unq_nm)))), 3)

  ratio_pep_pep <- round(rbind(Num_Peptide_Per_Pep,
                              Num_Cond_Peptide_Per_Pep,
                              Num_Rest_Peptide_Per_Pep,
                              Num_Rest_Peptide_Per_Pep / Num_Cond_Peptide_Per_Pep,
                              -log10(as.numeric(names(unq_pep)))), 3)

  rownames(ratio_pep_grp) <- c("Num_Peptide_Per_Grp", "Num_Cond_Peptide_Per_Grp", "Num_Rest_Peptide_Per_Grp",
                               "Num_Rest_Peptide_Per_Grp / Num_Cond_Peptide_Per_Grp", "-logP")
  rownames(ratio_pep_nm) <- c("Num_Peptide_Per_NM", "Num_Cond_Peptide_Per_NM", "Num_Rest_Peptide_Per_NM",
                               "Num_Rest_Peptide_Per_NM / Num_Cond_Peptide_Per_NM", "-logP")
  rownames(ratio_pep_pep) <- c("Num_Peptide_Per_Pep", "Num_Cond_Peptide_Per_Pep", "Num_Rest_Peptide_Per_Pep",
                              "Num_Rest_Peptide_Per_Pep / Num_Cond_Peptide_Per_Pep", "-logP")

  return(list(ans, ratio_pep_grp, ratio_pep_nm, ratio_pep_pep))
}
