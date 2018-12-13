#' Get background data frame (fasta library from Refseq).
#'
#' @param species A string for that the alignment is based on
#' which species, the options are 'human' and 'mouse'.
#'
#' @author Dongdong Zhan and Mengsha Tong
#'
#' @return A data frame of background
#' @export
#'
#' @examples
#' \dontrun{
#' background_df = get_global_background_df(species)
#' }


get_global_background_df <- function(species){
  requireNamespace('utils')

  # *** background ***
  cat('Reading background.\n')

  # Read fasta file: gi_fasta
  background_file_name <-  paste('STY_background_of_', species, '_for_motif_enrichment.txt', sep = '')
  BACKGROUND_FILE_PATH <- normalizePath(
    system.file(
      'extdata',
      'motif_library', species, background_file_name,
      package = "PhosMap"
    ),
    mustWork = F
  )
  if(!file.exists(BACKGROUND_FILE_PATH)){
    cat(BACKGROUND_FILE_PATH, ' -> ', 'No the file')
    stop('')
  }

  cat('Read background file of ', species, '.\n', sep = '')
  background_df <- utils::read.table(BACKGROUND_FILE_PATH, sep = '\t', header = T)
  cat('Read OK! ^_^', '\n')

  return(background_df)

}
