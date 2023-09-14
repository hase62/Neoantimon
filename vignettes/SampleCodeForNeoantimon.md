---
title: "SampleCodeForNeoantimon"
author: "T. Hasegawa"
date: "2020/6/27"
output:
pdf_document: default
html_document: default
---

<!--
%\VignetteEngine{knitr::knitr}
%\VignetteIndexEntry{Sample Code to Use Neoantimon}
-->

## Data Preparation and Sample Codes for Analysis

```r
#install.packages('devtools');
library(devtools);
```

```
## Loading required package: usethis
```

```r
install_github('hase62/Neoantimon');
```

```
## Downloading GitHub repo hase62/Neoantimon@HEAD
```

```
## askpass     (1.1   -> 1.2.0) [CRAN]
## dplyr       (1.1.2 -> 1.1.3) [CRAN]
## knitr       (1.43  -> 1.44 ) [CRAN]
## credentials (1.3.2 -> 2.0.1) [CRAN]
```

```
## Installing 4 packages: askpass, dplyr, knitr, credentials
```

```
## 
## The downloaded binary packages are in
## 	/var/folders/ns/yw46qhdj5cl8x7kx3z_mtjt00000gn/T//RtmpJZubD5/downloaded_packages
## ── R CMD build ──────────────────────────────────────────────────────────
##      checking for file ‘/private/var/folders/ns/yw46qhdj5cl8x7kx3z_mtjt00000gn/T/RtmpJZubD5/remotes69fb6786eedc/hase62-Neoantimon-5828ff3/DESCRIPTION’ ...  ✔  checking for file ‘/private/var/folders/ns/yw46qhdj5cl8x7kx3z_mtjt00000gn/T/RtmpJZubD5/remotes69fb6786eedc/hase62-Neoantimon-5828ff3/DESCRIPTION’
##   ─  preparing ‘Neoantimon’:
##      checking DESCRIPTION meta-information ...  ✔  checking DESCRIPTION meta-information
##   ─  checking for LF line-endings in source and make files and shell scripts
##   ─  checking for empty or unneeded directories
##   ─  looking to see if a ‘data/datalist’ file should be added
##   ─  building ‘Neoantimon_2.1.1.tar.gz’
##      
## 
```

```r
library(Neoantimon);
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("biomaRt")
```

```
## Bioconductor version 3.17 (BiocManager 1.30.22), R 4.3.1 (2023-06-16)
```

```
## Warning: package(s) not installed when version(s) same as or greater than current; use
##   `force = TRUE` to re-install: 'biomaRt'
```

```
## Old packages: 'KernSmooth', 'Matrix', 'RcppArmadillo', 'credentials',
##   'foreign', 'knitr', 'mgcv', 'minqa', 'nlme', 'reticulate', 'spatial',
##   'survival'
```

```r
library(biomaRt)
```
