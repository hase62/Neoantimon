## ----Preparation---------------------------------------------------------
print("install.packages('devtools')");
library(devtools);
print("install_github('hase62/Neoantimon')");
library(Neoantimon);

## ----Sample VCF file-----------------------------------------------------
data("sample_vcf")
print(sample_vcf, row.names = FALSE)

## ----Sample VCF (SV Fusion BND) file-------------------------------------
data("sample_sv_bnd")
print(sample_sv_bnd, row.names = FALSE)

## ----HLA Table for Class 1-----------------------------------------------
library(Neoantimon)
data("sample_hla_table_c1")
print(sample_hla_table_c1, row.names = FALSE)

## ----HLA Table for Class2------------------------------------------------
data("sample_hla_table_c2")
print(sample_hla_table_c2, row.names = FALSE)

## ----RNA Expression------------------------------------------------------
data("sample_rna_exp")
print(sample_rna_exp, row.names = FALSE)

## ----CopyNumber Information----------------------------------------------
data("sample_copynum")
print(sample_copynum, row.names = FALSE)

