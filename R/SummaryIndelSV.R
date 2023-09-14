#'Export Summary Count from Indel/SV Results
#'
#'@param Input Input file generated from MainSNVClass1,2.
#'
#'@param Mut_EL_th The threshold for mutant peptide to be neoantigen by EL.
#'
#'@param Mut_Rank_th The threshold for mutant peptide to be neoantigen by Rank.
#'
#'@param Total_RNA_th The total RNA expression threshold.
#'
#'@param Tumor_RNA_th The tumor specific RNA expression threshold.
#'
#'@param MutRatio_th The mutation ratio threshold.
#'
#'@param Weight The weight for alterations.
#'
#'@param WriteLongIndel If setting a file name, Write Long Indels of which the p-value is less than 0.05.
#'
#'@param IgnoreLongIndel Ignore Indels of which p-value is less than the indicated value for counting.
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
Export_Summary_IndelSV <- function(Input,
                                   Mut_EL_th = NA,
                                   Mut_Rank_th = NA,
                                   Total_RNA_th = NA,
                                   Tumor_RNA_th = NA,
                                   MutRatio_th = NA,
                                   Weight = NA,
                                   WriteLongIndel = NA,
                                   IgnoreLongIndel = 0,
                                   DupCount = FALSE){

  if((!is.na(Mut_EL_th) & !is.na(Mut_Rank_th)) |  (is.na(Mut_EL_th) & is.na(Mut_Rank_th))){
    print("Please Specify Either One of Mut_EL_th or Mut_Rank_th")
    return(NULL)
  }

  # EL or Rank
  m_th <- Mut_EL_th
  m_th_column <- "Mut_EL"
  if(!is.na(Mut_Rank_th)) {
    m_th <- Mut_Rank_th
    m_th_column <- "Mut_Rank"
  }

  # Attach Weight
  if(!is.na(Weight[1])) {
    Input <- cbind(Weight, Input)
    colnames(Input)[1] <- "Weight"
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
    Input[, match("Mutation_Position", index)] <- paste(Input[, match("HLA", index)], Input[, match("Mutation_Position", index)], sep = "_")
    Input[, match("Evaluated_Mutant_Peptide", index)] <- paste(Input[, match("HLA", index)], Input[, match("Evaluated_Mutant_Peptide", index)], sep = "_")
    num_of_hla <- length(unique(Input[, match("HLA", index)]))
  }

  # Count All
  Num_Alteration <- length(unique(Input[, match("Mutation_Position", index)]))
  Num_Peptide <- length(unique(Input[, match("Evaluated_Mutant_Peptide", index)]))

  # Conditioning
  if(IgnoreLongIndel > 0){
    Input <- Input[pvalues > IgnoreLongIndel, ]
  }
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
  Num_Cond_Alteration <- length(unique(Input[, match("Mutation_Position", index)]))
  Num_Cond_Peptide <- length(unique(Input[, match("Evaluated_Mutant_Peptide", index)]))

  # Extract by EL
  Input <- Input[as.numeric(Input[, match(m_th_column, index)]) < m_th, ]

  if(is.null(Input) | is.null(dim(Input)[1])) {
    if(length(Input) > 10) {
      Input <- t(Input)
    }
  }

  # Count Rest
  Num_Rest_Alteration <- length(unique(Input[, match("Mutation_Position", index)]))
  Num_Rest_Peptide <- length(unique(Input[, match("Evaluated_Mutant_Peptide", index)]))

  ans <- c(Num_Alteration / num_of_hla, Num_Cond_Alteration / num_of_hla, Num_Rest_Alteration / num_of_hla,
           Num_Peptide / num_of_hla, Num_Cond_Peptide / num_of_hla, Num_Rest_Peptide / num_of_hla)
  names(ans) <- c("Num_All_Alteration", "Num_Evaluated_Alteration", "Num_Alteration_Generating_NeoAg",
                  "Num_All_Peptide", "Num_Evaluated_Peptide", "Num_Peptide_Generating_NeoAg")

  if(!is.na(Weight[1])){
    alt_count <- sum(as.numeric(Input[match(unique(Input[, match("Mutation_Position", index)]),
                                            Input[, match("Mutation_Position", index)]),
                                      match("Weight", index)]))
    pep_count <- sum(as.numeric(Input[match(unique(Input[, match("Evaluated_Mutant_Peptide", index)]),
                                            Input[, match("Evaluated_Mutant_Peptide", index)]),
                                      match("Weight", index)]))
    ans <- c(ans[1:3], alt_count / num_of_hla, ans[4:6], pep_count / num_of_hla)
    names(ans)[c(4, 8)] <- c("Weighted_Num_Alteration_Generating_NeoAg", "Weighted_Num_Peptide_Generating_NeoAg")
  }
  return(ans)
}
