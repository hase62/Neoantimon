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

#' A Sample file for refFlat
#'
#' A dataset containing a part of refFlat data.
#'
#' @docType data
#' @keywords datasets
#' @name sample_refFlat.grch37
#' @usage data(sample_refFlat.grch37)
#' @format A data frame with 11 column.
NULL

#' A Sample file for refSeq RNA
#'
#' A dataset containing a part of refSeq RNA.
#'
#' @docType data
#' @keywords datasets
#' @name sample_refMrna.grch37.fa
#' @usage data(sample_refMrna.grch37.fa)
#' @format A data frame with 1 column.
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

#' A Format / Sample file for Annotated vcf file basef on Annovar.
#'
#' A dataset containing the variant information of a patient.
#'
#' @docType data
#' @keywords datasets
#' @name sample_vcf.annovar
#' @usage data(sample_vcf.annovar)
#' @format A data frame with 9 rows and variables including "Chr"	"Start"	"End"	"Ref"	"Alt"	"Func.refGene (exonic, intron, intergenic, ...)"	"ExonicFunc.refGene (exonic nonsynonymous, synonymous, insertion, ...)"	"AAChange.refGene (e.g., SLCO1C1:NM_001145944:exon7:c.692_693insG:p.L231fs ...)"
NULL

#' A Format / Sample file for Annotated vcf file based on VEP.
#'
#' A dataset containing the variant information of a patient.
#'
#' @docType data
#' @keywords datasets
#' @name sample_vcf.vep
#' @usage data(sample_vcf.vep)
#' @format A data frame with variables including "#Uploaded_variation"	"Location"	"Allele"	"Gene"	"Feature"	"Feature_type"	"Consequence"	"cDNA_position"	"CDS_position"	"Protein_position"	"Amino_acids	Codons"	"Existing_variation"	"Extra"
NULL

#' A Format / Sample file for snp informatin.
#'
#' A dataset containing snps information of a patient.
#'
#' @docType data
#' @keywords datasets
#' @name sample.snps.vcf
#' @usage data(sample.snps.vcf)
#' @format A data frame with variables including "#CHROM"	"POS"	"ID"	"REF"	"ALT"	"QUAL"	"FILTER"	"INFO"	"FORMAT".
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

#' Analyzed Result for SV Fusion CLASS1
#'
#' @docType data
#' @keywords result
#' @name sample_result_SVFusion_CLASS1_ALL
#' @usage data(sample_result_SVFusion_CLASS1_ALL)
NULL

#' Analyzed Result for SVFusion CLASS2
#'
#' @docType data
#' @keywords result
#' @name sample_result_SVFusion_CLASS2_ALL
#' @usage data(sample_result_SVFusion_CLASS2_ALL)
NULL

#' Analyzed Result for A DNA Fragment CLASS1
#'
#' @docType data
#' @keywords result
#' @name sample_result_SeqFragment_CLASS1_ALL
#' @usage data(sample_result_SeqFragment_CLASS1_ALL)
NULL

#' Analyzed Result for A DNA Fragment CLASS2
#'
#' @docType data
#' @keywords result
#' @name sample_result_SeqFragment_CLASS2_ALL
#' @usage data(sample_result_SeqFragment_CLASS2_ALL)
NULL
