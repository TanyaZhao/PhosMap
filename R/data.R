#' BRAFi data on quantification proteomics
#'
#' Data came from Ressa et al. MS experiments, they performed (phospho)proteomic profiling
#' of the WiDr colorectal cancer cells harboring the BRAF(V600E) mutation after treatment
#' using vemurafenib (BRAF inhibitor, abbr. BRAFi) in a time course at 0, 2, 6, 24, and 48 hour,
#' respectively.
#'
#' @name BRAFi
#'
#' @docType data
#'
#' @usage data(BRAFi)
#'
#' @format An list containning intermediate result from demo.
#'   \describe{
#'       \item{background_df}{A data frame for motif enrichment analysis as background.}
#'       \item{data_frame_normalization_with_control_no_pair}{A data frame containning phosphoproteomics data normalized by proteomics data.}
#'       \item{foreground_df}{A data frame for motif enrichment analysis as foreground.}
#'       \item{fuzzy_input_df}{A data frame for time course analysis as input.}
#'       \item{merge_df_with_phospho_peptides}{A merged phosphoproteomics data frame based on peptides files (unique ID).}
#'       \item{motif_group_m_ratio_df_mat}{A matrix for motif profile.}
#'       \item{phospho_data_normalization_and_filtering_STY}{A phosphoproteomics data frame after normalization and filtering.}
#'       \item{profiling_data_normalized}{A proteomics data frame after normalization and filtering.}
#'       \item{summary_df_of_unique_proteins_with_sites}{A data frame that phosphorylation sites had been mapping to protein sequence and eliminated redundancy.}
#'       \item{group}{A factor for experiment group information.}
#'  }
#'
#' @keywords datasets
#'
#' @references Ressa, A, et al. (2018) A System-wide Approach to Monitor Responses to Synergistic BRAF and EGFR Inhibition in Colorectal Cancer Cells, Molecular & cellular proteomics : MCP, 17, 1892-1908.
#'
#' @source \url{https://www.ebi.ac.uk/pride/archive/projects/PXD007740/}
#'
'BRAFi'




