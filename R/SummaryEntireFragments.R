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
Export_Summary_Entire_Fragments <- function(Input,
                                          Mut_EL_th = NA,
                                          Mut_Rank_th = NA,
                                          Total_RNA_th = NA,
                                          Tumor_RNA_th = NA,
                                          MutRatio_th = NA,
                                          WriteLongIndel = NA,
                                          DupCount = FALSE){

  if((!is.na(Mut_EL_th) & !is.na(Mut_Rank_th)) | (is.na(Mut_EL_th) & is.na(Mut_Rank_th))){
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

  # Calculate Pvalues
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

  # Duplication
  num_of_hla <- 1
  if(DupCount){
    Input[, match("Evaluated_Mutant_Peptide", index)] <- paste(Input[, match("HLA", index)],
                                                               Input[, match("Evaluated_Mutant_Peptide", index)], sep = "_")
    num_of_hla <- length(unique(Input[, match("HLA", index)]))
  }

  # Count All
  res <- NULL
  for(tag in c("GroupID", "NM_ID", "Mutant_Peptide")){
    Input_tag <- Input[, match(tag, index)]
    if(tag == "Mutant_Peptide"){
      Input_tag <- paste(Input_tag, Input[, match("Gene", index)], Input[, match("GroupID", index)], sep = "-")
    }
    unq_tag <- unique(Input_tag)
    id_unq_tag <- sapply(unq_tag, function(x) which(!is.na(match(Input_tag, x))))
    Input_unq_tag <- lapply(id_unq_tag, function(x) Input[x, ])
    names(unq_tag) <- lapply(id_unq_tag, function(x) min(as.numeric(pvalues[x])))

    eval_col <- match("Evaluated_Mutant_Peptide", index)
    Num_Peptide_Per_tag <- unlist(lapply(Input_unq_tag, function(x) length(unique(x[, eval_col]))))

    # Conditioning
    if(!is.na(Total_RNA_th)){
      rna_col <- match("Total_RNA", index)
      Input_unq_tag <- lapply(Input_unq_tag, function(x) {x[as.numeric(x[, rna_col]) > Total_RNA_th, ]
                                                          if(!is.matrix(x)) {t(x)}else{x}})
    }
    if(!is.na(Tumor_RNA_th)){
      rna_col <- match("Tumor_RNA", index)
      Input_unq_tag <- lapply(Input_unq_tag, function(x) {x[as.numeric(x[, rna_col]) > Tumor_RNA_th, ]
                                                          if(!is.matrix(x)) {t(x)}else{x}})
    }
    if(!is.na(MutRatio_th)){
      rna_col <- match("MutRatio", index)
      Input_unq_tag <- lapply(Input_unq_tag, function(x) {x[as.numeric(x[, rna_col]) > MutRatio_th, ]
                                                          if(!is.matrix(x)) {t(x)}else{x}})
    }

    # Count Conditioned
    Num_Cond_Peptide_Per_tag <- unlist(lapply(Input_unq_tag, function(x) length(unique(x[, eval_col]))))

    # Extract by EL
    m_col <- match(m_th_column, index)
    Input_unq_tag <- lapply(Input_unq_tag, function(x) {x <- x[as.numeric(x[, m_col]) < m_th, ]
                                                        if(!is.matrix(x)) {t(x)}else{x}})

    # Count Rest
    Num_Rest_Peptide_Per_tag <- unlist(lapply(Input_unq_tag, function(x) length(unique(x[, eval_col]))))

    ratio_pep_tag <- round(rbind(Num_Peptide_Per_tag / num_of_hla,
                                 Num_Cond_Peptide_Per_tag / num_of_hla,
                                 Num_Rest_Peptide_Per_tag / num_of_hla,
                                 Num_Rest_Peptide_Per_tag / Num_Cond_Peptide_Per_tag,
                                 -log10(as.numeric(names(unq_tag)))), 3)

    if(tag == "GroupID") {
      Input_Gene <- Input[, match("Gene", index)]
      unq_grp_gene <- Input_Gene[match(unq_tag, Input[, match("GroupID", index)])]
      unq_grp_gene <- sapply(unq_grp_gene, function(x) strsplit(x, "_")[[1]][2])
      colnames(ratio_pep_tag) <- unq_grp_gene
      rownames(ratio_pep_tag) <- c("Num_Peptide_Per_Grp",
                                   "Num_Cond_Peptide_Per_Grp",
                                   "Num_Rest_Peptide_Per_Grp",
                                   "Num_Rest_Peptide_Per_Grp / Num_Cond_Peptide_Per_Grp",
                                   "-logP")
    }

    if(tag == "NM_ID") {
      colnames(ratio_pep_tag) <- unq_tag
      rownames(ratio_pep_tag) <- c("Num_Peptide_Per_NM",
                                  "Num_Cond_Peptide_Per_NM", "Num_Rest_Peptide_Per_NM",
                                  "Num_Rest_Peptide_Per_NM / Num_Cond_Peptide_Per_NM",
                                  "-logP")
    }
    if(tag == "Mutant_Peptide") {
      colnames(ratio_pep_tag) <- unq_tag
      rownames(ratio_pep_tag) <- c("Num_Peptide_Per_Pep",
                                   "Num_Cond_Peptide_Per_Pep",
                                   "Num_Rest_Peptide_Per_Pep",
                                   "Num_Rest_Peptide_Per_Pep / Num_Cond_Peptide_Per_Pep",
                                   "-logP")
    }
    res <- c(res, list(ratio_pep_tag))
  }

  return(res)
}
