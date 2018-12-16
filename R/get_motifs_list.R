#' Motif enrichment using rmotifx.
#'
#' @param foreground A vector for aligned sequences as the foreground input.
#' @param background A vector for aligned sequences as the background input.
#' @param center_vector A vector for aligned centers.
#' @param motifx_pvalue A numeric value for selecting motifs that meets the minimum cutoff.
#'
#' @author Dongdong Zhan and Mengsha Tong
#'
#' @references Omar Wagih (2014). rmotifx: An iterative statistical approach to the discovery of biological sequence motifs. R package version 1.0.
#'
#' @return A list from results of motif enrichment.
#' @export
#'
#' @examples
#' \dontrun{
#' motifs_list = get_motifs_list(foreground, background, center_vector, motifx_pvalue)
#' }

get_motifs_list <- function(foreground, background, center_vector, motifx_pvalue){
  motifs_list <- list()
  center_vector_len <- length(center_vector)
  for(i in 1:center_vector_len){
    center <- center_vector[i]
    motifs <- get_motif_analysis_summary_list(foreground, background, center = center, min_sequence_count = 1, min_pvalue = motifx_pvalue)
    motifs_list[[i]] <- motifs
  }
  names(motifs_list) <- center_vector
  return(motifs_list)
}
