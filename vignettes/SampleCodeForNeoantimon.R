## ----Preparation--------------------------------------------------------------
#install.packages('devtools');
library(devtools);
install_github('hase62/Neoantimon');
library(Neoantimon);
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("biomaRt")
library(biomaRt)

