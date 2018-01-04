#' Neoantimon
#'
#' Calculate Lists of Candidate Neoantingens from SNVs, Indels, and SV fusions to MHC Class1 and Class2.
#' First use MainSNVClass1, MainSNVClass2, MainINDELClass1, MainINDELClass2, MainSVFUSIONClass1, and MainSVFUSIONClass2.
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

#' A Format / Sample file for Annotated vcf file.
#'
#' A dataset containing the variant information of a patient.
#'
#' @docType data
#' @keywords datasets
#' @name sample_vcf
#' @usage data(sample_vcf)
#' @format A data frame with 9 rows and variables including "Chr"	"Start"	"End"	"Ref"	"Alt"	"Func.refGene (exonic, intron, intergenic, ...)"	"ExonicFunc.refGene (exonic nonsynonymous, synonymous, insertion, ...)"	"AAChange.refGene (e.g., SLCO1C1:NM_001145944:exon7:c.692_693insG:p.L231fs ...)"
NULL

#' A Format / Sample file for Annotated vcf file.
#'
#' A dataset containing the variant information of a patient.
#'
#' @docType data
#' @keywords datasets
#' @name sample_sv_bnd
#' @usage data(sample_sv_bnd)
#' @format A data frame with 9 rows and variables including "Chr"	"Start"	"End"	"Ref"	"Alt (BND format)"	"Func.refGene (exonic, intron, intergenic, ...)"	"ExonicFunc.refGene (exonic nonsynonymous, synonymous, insertion, ...)"	"mateID (e.g., SVMERGE1_1)"
NULL

#' Analyzed Result for SNV CLASS1
#'
#' @docType data
#' @keywords result
#' @name sample_result_SNV_CLASS1_ALL
#' @usage data(sample_result_SNV_CLASS1_ALL)
NULL

#' Analyzed Result for SNV CLASS2
#'
#' @docType data
#' @keywords result
#' @name sample_result_SNV_CLASS2_ALL
#' @usage data(sample_result_SNV_CLASS2_ALL)
NULL

#' Analyzed Result for INDEL CLASS1
#'
#' @docType data
#' @keywords result
#' @name sample_result_INDEL_CLASS1_ALL
#' @usage data(sample_result_INDEL_CLASS1_ALL)
NULL

#' Analyzed Result for INDEL CLASS2
#'
#' @docType data
#' @keywords result
#' @name sample_result_INDEL_CLASS2_ALL
#' @usage data(sample_result_INDEL_CLASS2_ALL)
NULL
