#' Get motifs and their corresponding aligned sequences form from foreground.
#'
#' @param motifs_list A list from motif enrichment analysis.
#' @param foreground A vector for aligned sequences.
#'
#' @author Dongdong Zhan and Mengsha Tong
#'
#' @references Hadley Wickham (2018). stringr: Simple, Consistent Wrappers for Common String Operations. R package version 1.3.0.\
#' https://CRAN.R-project.org/package=stringr.
#'
#' @return A list containing motifs and the corresponding sequences from foreground.
#' @export
#'
#' @examples
#' \dontrun{
#' mots_match_list = get_foreground_sequences_mapped_to_motifs(
#'   motifs_list,
#'   foreground
#' )
#' }


get_foreground_sequences_mapped_to_motifs <- function(motifs_list, foreground){
  requireNamespace('stats')
  cat('\n', 'Find sequences in foreground that are mapped to specific motif. \n')
  requireNamespace('stringr')
  motifs_list_len <- length(motifs_list)
  mots_match_list <- list()
  mots_match_index <- 0
  mots_patterns_valid <- NULL
  for(mots_index in seq_len(motifs_list_len)){
    mots <- motifs_list[[mots_index]]
    mots_patterns <- as.vector(mots$motif)
    mots_patterns_len <- length(mots_patterns)
    for(i in seq_len(mots_patterns_len)){
      mots_pattern <- mots_patterns[i]
      regex_pattern <- stringr::str_replace_all(mots_pattern, '\\.', '[A-Z_]{1}')
      match_results <- apply(data.frame(foreground), 1, function(x, regex_pattern){
        x <- stringr::str_extract(x, regex_pattern)
      }, regex_pattern = regex_pattern)
      match_results <- as.vector(stats::na.omit(match_results))
      if(length(match_results)>0){
        mots_match_index <- mots_match_index + 1
        mots_match_list[[mots_match_index]] <- match_results
        mots_patterns_valid <- union(mots_patterns_valid, mots_patterns)
      }

    }
  }
  names(mots_match_list) <- mots_patterns_valid
  return(mots_match_list)
}
