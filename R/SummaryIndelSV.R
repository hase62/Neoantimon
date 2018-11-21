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
#'@param Weight_rf2 The weight for alterations.
#'
#'@param Weight_rf3 The weight for alterations.
#'
#'@param WriteLongIndel If setting a file name, Write Long Indels.
#'
#'@return Num_Alteration The number of evaluated alterations.
#'
#'@return Num_Alteration_Generating_NeoAg The number of evaluated alterations that can generate neoantigen.
#'
#'@return Num_Peptide The number of evaluated peptifdes.
#'
#'@return Num_Peptide_Generating_NeoAg The number of evaluated peptides that can be neoantigen.
#'
#'
#'@export
Export_Summary_IndelSV <- function(Input,
                               mut_IC50_th,
                               Total_RNA_th = NA,
                               Tumor_RNA_th = NA,
                               MutRatio_th = NA,
                               Weight_rf2 = NA,
                               Weight_rf3 = NA,
                               WriteLongIndel = NA){

  index <- colnames(Input)
  if(!is.na(WriteLongIndel)){
    theoretical_pro <- 3 * 4^3 / 61^2
    tmp <- Input[sapply(Input[, match("Mutant_Peptide", index)], function(x) (1-theoretical_pro)^nchar(x)) < 0.01, ]
    tmp <- match(unique(tmp[, grep("Mutation_Position", index)]),
                 Input[, grep("Mutation_Position", index)])
    if(length(tmp) > 0){
      write.table(Input[tmp, (match("Chr", index)):ncol(Input)],
                paste(WriteLongIndel, ".txt", sep = ""),
                row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
    }
  }

  if(!is.na(Weight_rf2[1]) && !is.na(Weight_rf3[1])){
    pos <- match(unique(Input[, grep("NM_ID", index)]), Input[, grep("NM_ID", index)])
    pos <- apply2(gsub("-", "", Input[pos, c(11, 12)]), 1, function(x) nchar(x[1]) - nchar(x[2])) %% 3
    Weight = ifelse(pos == 1, Weight_rf2, ifelse(pos == 2, Weight_rf3, 0))
    Input <- cbind(Input, match(Input[, match("NM_ID", index)], unique(Input[, grep("NM_ID", index)])))
    print("Set Weight as")
    print(paste(unique(Input[, grep("NM_ID", index)]), "is", Weight))
  }
  Num_Alteration <-length(unique(Input[, grep("Mutation_Position", index)]))
  Num_Peptide <-length(unique(Input[, grep("Evaluated_Mutant", index)]))

  Input <- Input[as.numeric(Input[, match("Mut_IC50", index)]) < mut_IC50_th, ]
  if(!is.na(Total_RNA_th)){
    Input <- Input[as.numeric(Input[, match("Total_RNA", index)]) > Total_RNA_th, ]
  }
  if(!is.na(Tumor_RNA_th)){
    Input <- Input[as.numeric(Input[, match("Tumor_RNA", index)]) > Tumor_RNA_th, ]
  }
  if(!is.na(MutRatio_th)){
    Input <- Input[as.numeric(Input[, match("MutRatio", index)]) > MutRatio_th, ]
  }

  if(is.null(Input)) {
    if(length(Input) > 10) {
      Input <- t(Input)
    } else {
      return(NULL)
    }
  }
  if(is.na(Weight_rf2[1]) | is.na(Weight_rf3[1])){
    alt_count <- length(unique(Input[, grep("Mutation_Position", index)]))
    pep_count <- length(unique(Input[, grep("Evaluated_Mutant", index)]))
  } else {
    alt_count <- sum(Weight[as.numeric(Input[match(unique(Input[, grep("Mutation_Position", index)]),
                      Input[, grep("Mutation_Position", index)]), ncol(Input)])])
    pep_count <- sum(Weight[as.numeric(Input[match(unique(Input[, grep("Evaluated_Mutant", index)]),
                                                   Input[, grep("Evaluated_Mutant", index)]), ncol(Input)])])
  }

  ans <- c(Num_Alteration, alt_count, Num_Peptide, pep_count)
  names(ans) <- c("Num_Alteration", "Num_Alteration_Generating_NeoAg",
                  "Num_Peptide", "Num_Peptide_Generating_NeoAg")

  return(ans)
}
