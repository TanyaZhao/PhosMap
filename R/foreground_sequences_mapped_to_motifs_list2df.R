#' Convert the list that consists of motifs and the corresponding sequences to data frame.
#'
#' @param foreground_sequences_mapped_to_motifs A list that consists of motifs and the corresponding sequences.
#' @author Dongdong Zhan and Mengsha Tong
#'
#' @return A data frame that consist of aligned sequences and the corresponding motifs.
#' @export
#'
#' @examples
#' \dontrun{
#' df <- foreground_sequences_mapped_to_motifs_list2df(
#'   foreground_sequences_mapped_to_motifs
#' )
#' }

foreground_sequences_mapped_to_motifs_list2df <- function(
  foreground_sequences_mapped_to_motifs
){
  foreground_sequences_mapped_to_motifs_len <-  length(foreground_sequences_mapped_to_motifs)
  motif_rep_v <- NULL
  map_seq_v <- NULL
  for(i in seq_len(foreground_sequences_mapped_to_motifs_len)){
    motif <- names(foreground_sequences_mapped_to_motifs[i])
    map_seq <- foreground_sequences_mapped_to_motifs[[i]]
    motif_rep <- rep(motif, length(map_seq))
    motif_rep_v <- c(motif_rep_v, motif_rep)
    map_seq_v <- c(map_seq_v, map_seq)
  }
  df <- data.frame(map_seq_v, motif_rep_v)
  colnames(df) <- c('Aligned_seq', 'Motif')
  return(df)
}





