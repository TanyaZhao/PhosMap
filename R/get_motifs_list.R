#' Motif enrichment using rmotifx.
#'
#' @param foreground A vector for aligned sequences as the foreground input.
#' @param background A vector for aligned sequences as the background input.
#' @param central_vector A vector for aligned centers.
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
#' motifs_list = get_motifs_list(foreground, background, central_vector, motifx_pvalue)
#' }


get_motifs_list <- function(foreground, background, central_vector, motifx_pvalue){
  requireNamespace('rmotifx')
  motifs_list <- list()
  central_vector_len <- length(central_vector)
  for(i in 1:central_vector_len){
    central <- central_vector[i]
    motifs <- rmotifx::motifx(foreground, background, central.res = central, min.seqs = 1, pval.cutoff = 1e-2)
    motifs_list[[i]] <- motifs
  }
  names(motifs_list) <- central_vector
  return(motifs_list)
}
