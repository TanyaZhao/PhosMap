% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_all_psites_and_normalize_data_with_FOT5.R
\name{get_all_psites_and_normalize_data_with_FOT5}
\alias{get_all_psites_and_normalize_data_with_FOT5}
\title{To normalize data and filter data only including phosphorylation sites.}
\usage{
get_all_psites_and_normalize_data_with_FOT5(data_frame,
  experiment_code_file_path)
}
\arguments{
\item{data_frame}{A data frame containing IDs and quantification values merged from multi-experiments as input.}

\item{experiment_code_file_path}{A file path of storing experiment codes as input. The experiment codes are required to keep pace with column names of Value.}
}
\value{
A data frame after normalization and filtering (x 1e5)
}
\description{
To normalize data and filter data only including phosphorylation sites.
}
\examples{
\dontrun{
ptypes_df = get_normalized_data_FOT5(data_frame, experiment_code_file_path)
}
}
