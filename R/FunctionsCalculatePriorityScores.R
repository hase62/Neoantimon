#'Calculate priority scores
#'
#'@param result A results file generated from main functions
#
#'@param useRNAvaf To indicate whether this function uses DNA VAF or RNA VAF for the calculation.
#'
#'@return P_I: Priority score using the EL.
#'
#'@return P_R: Priority score using the percentage of rank affinity.
#'
#'@return P: Priority score implemented in MuPeXI (Bjerregaard et al. 2017).
#'
#'@export
CalculatePriorityScores <- function(result, useRNAvaf = FALSE){
  if(useRNAvaf){
    vaf <- result[, match("Tumor_RNA_Ratio", colnames(result))]
    vaf <- sapply(vaf, function(x) as.numeric(strsplit(x, "/")[[1]][1]) / as.numeric(strsplit(x, "/")[[1]][2]))
  } else {
    vaf <- as.numeric(result[, match("Tumor_Depth", colnames(result))]) / as.numeric(result[, match("Total_Depth", colnames(result))])
  }
  vaf <- ifelse(is.na(vaf), 1.0, vaf)
  gA <- vaf * tanh(5 * vaf)
  E <- as.numeric(result[, match("Total_RNA", colnames(result))])

  j1 <- match("Evaluated_Mutant_Peptide", colnames(result))
  j2 <- match("Evaluated_Wt_Peptide", colnames(result))
  if(!is.na(j2)) {
    M <- sapply(1:nrow(result),
              function(i){
                m <- strsplit(result[i, j1], "")[[1]]
                w <- strsplit(result[i, j2], "")[[1]]
                if(length(w) - length(m) < 0){
                  M <- 2^-length(m)
                } else {
                  M <- 2^-min(sapply(seq(1, length(w) - length(m) + 1, 1),
                                     function(point) length(m) - length(which(m == w[point:(point + length(m) - 1)]))))
                }
              return(M)
              })
  } else {
    M <- 0
  }

  L_I_M <- 1 / (1 + exp(0.015 * (as.numeric(result[, match("Mut_EL", colnames(result))]) - 500)))
  L_R_M <- 1 / (1 + exp(5 * (as.numeric(result[, match("Mut_Rank", colnames(result))]) - 2)))
  j_3 <- match("Wt_EL", colnames(result))
  if(!is.na(j_3)){
    L_I_W <- 1 / (1 + exp(0.015 * (as.numeric(result[, match("Wt_EL", colnames(result))]) - 500)))
    L_R_W <- 1 / (1 + exp(5 * (as.numeric(result[, match("Wt_Rank", colnames(result))]) - 2)))
  } else {
    L_I_W <- 0
    L_R_W <- 0
  }

  j_4 <- match("MutRatio", colnames(result))
  C <- ifelse(!is.na(result[, j_4]) & result[, j_4] != "NA", 1 / (1 - exp(-30 * (as.numeric(result[, j_4]) - 0.8))), rep(1, nrow(result)))

  P_I <-  L_I_M * tanh(gA * E) * (1 - M * L_I_W) * C
  P_R <-  L_R_M * tanh(gA * E) * (1 - M * L_R_W) * C
  P <- L_R_M * vaf * tanh(E) * (1 - M * L_R_W)

  pri_scores <- cbind(P_I, P_R, P)
  colnames(pri_scores) <- c("P_I", "P_R", "P")
  rownames(pri_scores) <- NULL

  return(cbind(result, pri_scores))
}
