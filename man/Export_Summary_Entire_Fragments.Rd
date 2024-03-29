% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SummaryEntireFragments.R
\name{Export_Summary_Entire_Fragments}
\alias{Export_Summary_Entire_Fragments}
\title{Export Summary Count from Indel/SV Results}
\usage{
Export_Summary_Entire_Fragments(
  Input,
  Mut_EL_th = NA,
  Mut_Rank_th = NA,
  Total_RNA_th = NA,
  Tumor_RNA_th = NA,
  MutRatio_th = NA,
  WriteLongIndel = NA,
  DupCount = FALSE
)
}
\arguments{
\item{Input}{Input file generated from MainSNVClass1,2.}

\item{Mut_EL_th}{The threshold for mutant peptide to be neoantigen by EL.}

\item{Mut_Rank_th}{The threshold for mutant peptide to be neoantigen by Rank.}

\item{Total_RNA_th}{The total RNA expression threshold.}

\item{Tumor_RNA_th}{The tumor specific RNA expression threshold.}

\item{MutRatio_th}{The mutation ratio threshold.}

\item{WriteLongIndel}{If setting a file name, Write Long Indels of which the p-value is less than 0.05.}

\item{DupCount}{Count for each different HLA type}
}
\value{
Num_Alteration The number of evaluated alterations.

Num_Alteration_Generating_NeoAg The number of evaluated alterations that can generate neoantigen.

Num_Peptide The number of evaluated peptifdes.

Num_Peptide_Generating_NeoAg The number of evaluated peptides that can be neoantigen.
}
\description{
Export Summary Count from Indel/SV Results
}
