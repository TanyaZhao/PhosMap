% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_motifs_list.R
\name{get_motifs_list}
\alias{get_motifs_list}
\title{Motif enrichment using rmotifx.}
\usage{
get_motifs_list(foreground, background, center_vector, motifx_pvalue)
}
\arguments{
\item{foreground}{A vector for aligned sequences as the foreground input.}

\item{background}{A vector for aligned sequences as the background input.}

\item{center_vector}{A vector for aligned centers.}

\item{motifx_pvalue}{A numeric value for selecting motifs that meets the minimum cutoff.}
}
\value{
A list from results of motif enrichment.
}
\description{
Motif enrichment using rmotifx.
}
\examples{
\dontrun{
motifs_list = get_motifs_list(foreground, background, center_vector, motifx_pvalue)
}
}
\references{
Omar Wagih (2014). rmotifx: An iterative statistical approach to the discovery of biological sequence motifs. R package version 1.0.
}
\author{
Dongdong Zhan and Mengsha Tong
}
