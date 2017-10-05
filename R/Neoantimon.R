#' Neoantimon
#'
#' Calculate Lists of Candidate Neoantingens on SNVs and Indels to MHC Class1 and Class2.
#' First use MainSNVClass1, MainSNVClass2, MainINDELClass1, and MainINDELClass2.
#'
#' @name Neoantimon
#' @docType package
NULL

#' A Format / Sample file for HLA CLASS1 Table
#'
#' A dataset containing the HLA types of patients in each row.
#'
#' @docType data
#' @keywords datasets
#' @name sample_hla_table_c1
#' @usage data(sample_hla_table_c1)
#' @format A data frame with 3 rows and at most 7 variables
NULL

#' A Format / Sample file for HLA CLASS2 Table
#'
#' A dataset containing the HLA types of patients in each row.
#'
#' @docType data
#' @keywords datasets
#' @name sample_hla_table_c2
#' @usage data(sample_hla_table_c2)
#' @format A data frame with at least 3 row and at most 10 variables
NULL

#' A Format / Sample file for Copy Number Information
#'
#' A dataset containing the copy number information obtained by, e.g., ASCAT.
#'
#' @docType data
#' @keywords datasets
#' @name sample_copynum
#' @usage data(sample_copynum)
#' @format A data frame with 7 rows and 9 variables
NULL

#' A Format / Sample file for RNA Expression Information
#'
#' A dataset containing the RNA expression amount of patient for each gene.
#'
#' @docType data
#' @keywords datasets
#' @name sample_rna_exp
#' @usage data(sample_rna_exp)
#' @format A data frame with 22 rows and 3 variables
NULL

#' A Format / Sample file for Analyzed vcf file.
#'
#' A dataset containing the variant information of a patient.
#'
#' @docType data
#' @keywords datasets
#' @name sample_vcf
#' @usage data(sample_vcf)
#' @format A data frame with 9 rows and 38 variables
NULL
