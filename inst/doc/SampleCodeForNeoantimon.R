## ----Preparation---------------------------------------------------------
#install.packages('devtools');
library(devtools);
install_github('hase62/Neoantimon');
library(Neoantimon);

## ----Sample VCF file-----------------------------------------------------
data("sample_vcf")
head(sample_vcf, row.names = FALSE)

## ----Sample VCF (SV Fusion BND) file-------------------------------------
data("sample_sv_bnd")
head(sample_sv_bnd, row.names = FALSE)

## ----HLA Table for Class 1-----------------------------------------------
library(Neoantimon)
data("sample_hla_table_c1")
head(sample_hla_table_c1, row.names = FALSE)

## ----HLA Table for Class2------------------------------------------------
data("sample_hla_table_c2")
head(sample_hla_table_c2, row.names = FALSE)

## ----RNA Expression------------------------------------------------------
data("sample_rna_exp")
head(sample_rna_exp, row.names = FALSE)

## ----CopyNumber Information----------------------------------------------
data("sample_copynum")
head(sample_copynum, row.names = FALSE)

## ----Sample.Result.SNV.HLACLASS1-----------------------------------------
data("sample_result_SNV_CLASS1_ALL")
head(sample_result_SNV_CLASS1_ALL, row.names = FALSE)

## ----Sample.Result.SNV.HLACLASS2-----------------------------------------
data("sample_result_SNV_CLASS2_ALL")
head(sample_result_SNV_CLASS2_ALL, row.names = FALSE)

## ----Sample.Result.INDEL.HLACLASS1---------------------------------------
data("sample_result_INDEL_CLASS1_ALL")
head(sample_result_INDEL_CLASS1_ALL, row.names = FALSE)

## ----Sample.Result.INDEL.HLACLASS2---------------------------------------
data("sample_result_INDEL_CLASS2_ALL")
head(sample_result_INDEL_CLASS2_ALL, row.names = FALSE)

## ----Get Sample and Test Analysis----------------------------------------
TestAnalysis

