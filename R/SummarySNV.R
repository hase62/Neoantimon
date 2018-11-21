#'Export Summary Count from SNV Results
#'
#'@param Input Input file generated from MainSNVClass1,2.
#'
#'@param mut_IC50_th The threshold for mutant peptide to be neoantigen.
#'
#'@param wt_IC50_th The threshold for wt peptide to be neoantigen.
#'
#'@param Total_RNA_th The total RNA expression threshold.
#'
#'@param Tumor_RNA_th The tumor specific RNA expression threshold.
#'
#'@param MutRatio_th The mutation ratio threshold.
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
                               mut_IC50_th,
                               wt_IC50_th = NA,
                               Total_RNA_th = NA,
                               Tumor_RNA_th = NA,
                               MutRatio_th = NA){

  index <- colnames(Input)
  Num_Alteration <-length(unique(Input[, grep("Mutation_Position", index)]))
  Num_Peptide <-length(unique(Input[, grep("Evaluated_Mutant", index)]))

  Input <- Input[as.numeric(Input[, match("Mut_IC50", index)]) < mut_IC50_th, ]
  if(!is.na(wt_IC50_th)){
    Input <- Input[as.numeric(Input[, match("Wt_IC50", index)]) < wt_IC50_th, ]
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

  if(is.null(Input) | is.null(dim(Input)[1])) {
    if(length(Input) > 10) {
      Input <- t(Input)
    } else {
      return(NULL)
    }
  }

  alt_count <- length(unique(Input[, grep("Mutation_Position", index)]))
  pep_count <- length(unique(Input[, grep("Evaluated_Mutant", index)]))

  ans <- c(Num_Alteration, alt_count, Num_Peptide, pep_count)
  names(ans) <- c("Num_Alteration", "Num_Alteration_Generating_NeoAg",
                  "Num_Peptide", "Num_Peptide_Generating_NeoAg")

  return(ans)
}
