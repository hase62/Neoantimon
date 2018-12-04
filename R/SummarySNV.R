#'Export Summary Count from SNV Results
#'
#'@param Input Input file generated from MainSNVClass1,2.
#'
#'@param Mut_IC50_th The threshold for mutant peptide to be neoantigen.
#'
#'@param Mut_Rank_th The threshold for mutant peptide to be neoantigen.
#'
#'@param Wt_IC50_th The threshold for wt peptide to be neoantigen.
#'
#'@param Wt_Rank_th The threshold for wt peptide to be neoantigen.
#'
#'@param Total_RNA_th The total RNA expression threshold.
#'
#'@param Tumor_RNA_th The tumor specific RNA expression threshold.
#'
#'@param MutRatio_th The mutation ratio threshold.
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
Export_Summary_SNV <- function(Input,
                               Mut_IC50_th = NA,
                               Mut_Rank_th = NA,
                               Wt_IC50_th = NA,
                               Wt_Rank_th = NA,
                               Total_RNA_th = NA,
                               Tumor_RNA_th = NA,
                               MutRatio_th = NA,
                               DupCount = FALSE){

  if((!is.na(Mut_IC50_th) & !is.na(Mut_Rank_th)) |  (is.na(Mut_IC50_th) & is.na(Mut_Rank_th))){
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

  index <- colnames(Input)

  # Duplication
  num_of_hla <- 1
  if(DupCount){
    Input[, match("Mutation_Position", index)] <- paste(Input[, match("HLA", index)], Input[, match("Mutation_Position", index)], sep = "_")
    Input[, match("Evaluated_Mutant_Peptide", index)] <- paste(Input[, match("HLA", index)], Input[, match("Evaluated_Mutant_Peptide", index)], sep = "_")
    num_of_hla <- length(unique(Input[, match("HLA", index)]))
  }

  # Count All
  Num_Alteration <-length(unique(Input[, match("Mutation_Position", index)]))
  Num_Peptide <-length(unique(Input[, match("Evaluated_Mutant_Peptide", index)]))

  # Conditioning
  if(!is.na(Wt_IC50_th)){
    Input <- Input[as.numeric(Input[, match("Wt_IC50", index)]) > Wt_IC50_th, ]
  }
  if(!is.na(Wt_Rank_th)){
    Input <- Input[as.numeric(Input[, match("Wt_Rank", index)]) > Wt_Rank_th, ]
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
  Num_Cond_Alteration <-length(unique(Input[, match("Mutation_Position", index)]))
  Num_Cond_Peptide <-length(unique(Input[, match("Evaluated_Mutant_Peptide", index)]))

  # Extract by IC50
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

  return(ans)
}
