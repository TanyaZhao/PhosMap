% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_summary_df_list_from_ksea.R
\name{get_summary_df_list_from_ksea}
\alias{get_summary_df_list_from_ksea}
\title{Get a data frame only containing inforamtion of kinase inferred by KSEA}
\usage{
get_summary_df_list_from_ksea(ptypes_data, species = "human",
  log2_label = TRUE, ratio_cutoff = 3)
}
\arguments{
\item{ptypes_data}{A data frame of phosphorylation data after normalization.}

\item{species}{A string representing the species of imported data, the options are human, mouse and rat.}

\item{log2_label}{A boolean value representing whether data is logarithmetics, the default is FALSE.}

\item{ratio_cutoff}{A cutoff that depicts quantification changes at phosphorylation level relative to profiling level, the default is 3.}
}
\value{
A data frame that consists of kinase, psite, substrate, counting byond ratio_cutoff and corresponding original value.
}
\description{
Get a data frame only containing inforamtion of kinase inferred by KSEA
}
\examples{
\dontrun{
summary_df_list_from_ksea = get_summary_df_list_from_ksea(
  ptypes_data,
  species = 'human',
  log2_label = TRUE,
  ratio_cutoff = 3
)
}



}
\author{
Dongdong Zhan and Mengsha Tong
}
