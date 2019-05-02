#'Export Summary Count from Indel/SV Results
#'
#'@param Input Input file generated from MainSNVClass1,2.
#'
#'@param Mut_IC50_th The threshold for mutant peptide to be neoantigen by IC50.
#'
#'@param Mut_Rank_th The threshold for mutant peptide to be neoantigen by Rank.
#'
#'@param Total_RNA_th The total RNA expression threshold.
#'
#'@param Tumor_RNA_th The tumor specific RNA expression threshold.
#'
#'@param MutRatio_th The mutation ratio threshold.
#'
#'@param WriteLongIndel If setting a file name, Write Long Indels of which the p-value is less than 0.05.
#'
#'@param DupCount Count for each different HLA type
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
                                     Mut_IC50_th = NA,
                                     Mut_Rank_th = NA,
                                     Total_RNA_th = NA,
                                     Tumor_RNA_th = NA,
                                     MutRatio_th = NA,
                                     WriteLongIndel = NA,
                                     DupCount = FALSE){
  
  if((!is.na(Mut_IC50_th) & !is.na(Mut_Rank_th)) | (is.na(Mut_IC50_th) & is.na(Mut_Rank_th))){
    print("Please Specify Either One of Mut_IC50_th or Mut_Rank_th")
    return(NULL)
  }
  
  # IC50 or Rank
  m_th <- Mut_IC50_th
  m_th_column <- "Mut_IC50"
  if(!is.na(Mut_Rank_th)) {
    m_th <- Mut_Rank_th
    m_th_column <- "Mut_Rank"
  }
  
  # Calculate Pvalues
  theoretical_pro <- 3 * 4^3 / 61^2
  pvalues <- sapply(Input[, match("Mutant_Peptide", colnames(Input))], function(x) (1 - theoretical_pro)^nchar(x))
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
  
  # Duplication
  num_of_hla <- 1
  if(DupCount){
    Input[, match("Evaluated_Mutant_Peptide", index)] <- paste(Input[, match("HLA", index)], Input[, match("Evaluated_Mutant_Peptide", index)], sep = "_")
    num_of_hla <- length(unique(Input[, match("HLA", index)]))
  }
  
  # Count All
  unq_grp <- unique(Input[, match("GroupID", index)])
  unq_grp_gene <- Input[match(unq_grp, Input[, match("GroupID", index)]), match("Gene", index)]
  unq_grp_gene <- sapply(unq_grp_gene, function(x) strsplit(x, "_")[[1]][2])
  unq_nm <- unique(Input[, match("NM_ID", index)])
  unq_pep <- unique(Input[, match("Mutant_Peptide", index)])
  names(unq_grp) <- sapply(unq_grp, function(x) min(as.numeric(Input[which(!is.na(match(Input[, match("GroupID", index)], x))), match("Pvalue", index)])))
  names(unq_nm) <- sapply(unq_nm, function(x) min(as.numeric(Input[which(!is.na(match(Input[, match("NM_ID", index)], x))), match("Pvalue", index)])))
  names(unq_pep) <- sapply(unq_pep, function(x) min(as.numeric(Input[which(!is.na(match(Input[, match("Mutant_Peptide", index)], x))), match("Pvalue", index)])))
  
  Num_Peptide_Per_Grp <- sapply(unq_grp, function(x) length(unique(Input[!is.na(match(Input[, match("GroupID", index)], x)), match("Evaluated_Mutant_Peptide", index)])))
  Num_Peptide_Per_NM <- sapply(unq_nm, function(x) length(unique(Input[!is.na(match(Input[, match("NM_ID", index)], x)), match("Evaluated_Mutant_Peptide", index)])))
  Num_Peptide_Per_Pep <- sapply(unq_pep, function(x) length(unique(Input[!is.na(match(Input[, match("Mutant_Peptide", index)], x)), match("Evaluated_Mutant_Peptide", index)])))
  
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
  Num_Cond_Peptide_Per_Grp <- sapply(unq_grp, function(x) length(unique(Input[!is.na(match(Input[, match("GroupID", index)], x)), match("Evaluated_Mutant_Peptide", index)])))
  Num_Cond_Peptide_Per_NM <- sapply(unq_nm, function(x) length(unique(Input[!is.na(match(Input[, match("NM_ID", index)], x)), match("Evaluated_Mutant_Peptide", index)])))
  Num_Cond_Peptide_Per_Pep <- sapply(unq_pep, function(x) length(unique(Input[!is.na(match(Input[, match("Mutant_Peptide", index)], x)), match("Evaluated_Mutant_Peptide", index)])))
  
  # Extract by IC50
  Input <- Input[as.numeric(Input[, match(m_th_column, index)]) < m_th, ]
  
  if(is.null(Input) | is.null(dim(Input)[1])) {
    if(length(Input) > 10) {
      Input <- t(Input)
    }
  }
  
  # Count Rest
  Num_Rest_Peptide_Per_Grp <- sapply(unq_grp, function(x) length(unique(Input[!is.na(match(Input[, match("GroupID", index)], x)), match("Evaluated_Mutant_Peptide", index)])))
  Num_Rest_Peptide_Per_NM <- sapply(unq_nm, function(x) length(unique(Input[!is.na(match(Input[, match("NM_ID", index)], x)), match("Evaluated_Mutant_Peptide", index)])))
  Num_Rest_Peptide_Per_Pep <- sapply(unq_pep, function(x) length(unique(Input[!is.na(match(Input[, match("Mutant_Peptide", index)], x)), match("Evaluated_Mutant_Peptide", index)])))
  
  ratio_pep_grp <- round(rbind(Num_Peptide_Per_Grp / num_of_hla,
                               Num_Cond_Peptide_Per_Grp / num_of_hla,
                               Num_Rest_Peptide_Per_Grp / num_of_hla,
                               Num_Rest_Peptide_Per_Grp / Num_Cond_Peptide_Per_Grp,
                               -log10(as.numeric(names(unq_grp)))), 3)
  
  ratio_pep_nm <- round(rbind(Num_Peptide_Per_NM / num_of_hla,
                              Num_Cond_Peptide_Per_NM / num_of_hla,
                              Num_Rest_Peptide_Per_NM / num_of_hla,
                              Num_Rest_Peptide_Per_NM / Num_Cond_Peptide_Per_NM,
                              -log10(as.numeric(names(unq_nm)))), 3)
  
  ratio_pep_pep <- round(rbind(Num_Peptide_Per_Pep / num_of_hla,
                               Num_Cond_Peptide_Per_Pep / num_of_hla,
                               Num_Rest_Peptide_Per_Pep / num_of_hla,
                               Num_Rest_Peptide_Per_Pep / Num_Cond_Peptide_Per_Pep,
                               -log10(as.numeric(names(unq_pep)))), 3)
  
  colnames(ratio_pep_grp) <- unq_grp_gene
  rownames(ratio_pep_grp) <- c("Num_Peptide_Per_Grp", "Num_Cond_Peptide_Per_Grp", "Num_Rest_Peptide_Per_Grp",
                               "Num_Rest_Peptide_Per_Grp / Num_Cond_Peptide_Per_Grp", "-logP")
  colnames(ratio_pep_nm) <- unq_nm
  rownames(ratio_pep_nm) <- c("Num_Peptide_Per_NM", "Num_Cond_Peptide_Per_NM", "Num_Rest_Peptide_Per_NM",
                              "Num_Rest_Peptide_Per_NM / Num_Cond_Peptide_Per_NM", "-logP")
  colnames(ratio_pep_pep) <- unq_pep
  rownames(ratio_pep_pep) <- c("Num_Peptide_Per_Pep", "Num_Cond_Peptide_Per_Pep", "Num_Rest_Peptide_Per_Pep",
                               "Num_Rest_Peptide_Per_Pep / Num_Cond_Peptide_Per_Pep", "-logP")
  
  return(list(ratio_pep_grp, ratio_pep_nm, ratio_pep_pep))
}
