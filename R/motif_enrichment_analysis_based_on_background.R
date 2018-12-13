#' Motif enrichment based on global background (fasta library from Refseq).
#'
#' @param foreground A vector for aligned sequence of foreground.
#' @param AA_in_protein A vector for the location of S/T/Y in sequence of protein.
#' @param background A vector for aligned sequence of background.
#' @param motifx_pvalue A numeric value for selecting motifs that meets the minimum cutoff.
#'
#' @author Dongdong Zhan and Mengsha Tong
#'
#' @return A list containing motifs and the corresponding sequences
#' @export
#'
#' @examples
#' \dontrun{
#' motifs_list = motif_enrichment_analysis_based_on_background(
#'   foreground,
#'   AA_in_protein,
#'   background,
#'   motifx_pvalue
#' )
#' }



# foreground = as.vector(foreground_df$aligned_seq)
# AA_in_protein = as.vector(foreground_df$AA_in_protein)
motif_enrichment_analysis_based_on_background <- function(foreground, AA_in_protein, background, motifx_pvalue){
  # foreground = as.vector(foreground)
  # background = as.vector(background$Aligned_Seq)
  central_vector_candidate <- c('S', 'T', 'Y')
  central_vector_candidate_len <- length(central_vector_candidate)
  central_vector <- NULL
  for(i in seq_len(central_vector_candidate_len)){
    central <- central_vector_candidate[i]
    if(length(grep(central, AA_in_protein)) > 0){
      central_vector <- c(central_vector, central)
    }
  }
  cat('Start executing motifx and find motif pattern. \n')
  cat('Foreground sequences: ', length(foreground), '.\n', sep = '')
  cat('Background sequences: ', length(background), '.\n', sep = '')
  cat('Phosphorylation: [', central_vector, '] exists in foreground.\n', sep = '')
  cat('Motifx pvalue cutoff: ', motifx_pvalue, '.\n', sep = '')
  motifs_list <- get_motifs_list(foreground, background, central_vector, motifx_pvalue)
  cat('Motifx analysis OK! ^_^', '\n')
  print(motifs_list)
  cat('\n')
  return(motifs_list)
}













































