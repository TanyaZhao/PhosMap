% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/check_foreground_and_background.R
\name{check_foreground_and_background_and_center}
\alias{check_foreground_and_background_and_center}
\title{Check input for motif analysis}
\usage{
check_foreground_and_background_and_center(foreground, background, center)
}
\arguments{
\item{foreground}{A vector for AA sequences with fixed length as foreground input.}

\item{background}{A vector for AA sequences with fixed length as background input.}

\item{center}{A character for center of k-mer.}
}
\value{
A list passing check steps
}
\description{
Check input for motif analysis
}
\examples{
\dontrun{
check_result_list <- check_foreground_and_background_and_center(
  foreground,
  background,
  center
)
}
}
\author{
Dongdong Zhan and Mengsha Tong
}
