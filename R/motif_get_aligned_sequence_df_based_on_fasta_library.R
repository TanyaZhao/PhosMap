#' Taking S/T/Y as the center, align sequence to fasta library by specific length.
#'
#' @param GI A vector for protein gi.
#' @param Sequence A vector for sequence of peptide.
#' @param AA_in_protein A vector for the location of S/T/Y in sequence of protein.
#' @param fixed_length A numeric value for aligned sequence,the default is 15.
#' @param species A string for that the alignment is based on which species, the options are 'human' and 'mouse'
#'
#' @author Dongdong Zhan and Mengsha Tong
#'
#' @references Hadley Wickham (2018). stringr: Simple, Consistent Wrappers for Common String Operations. R package version 1.3.0.\
#' https://CRAN.R-project.org/package=stringr.
#'
#' @return A data frame containing GI, Sequence, AA_in_protein, aligned_seq.
#' @export
#'
#' @examples
#' \dontrun{
#' aligned_sequence_df_based_on_fasta_library =
#'  motif_get_aligned_sequence_df_based_on_fasta_library(
#'    GI,
#'    Sequence,
#'    AA_in_protein,
#'    fixed_length,
#'    species
#'  )
#' }

motif_get_aligned_sequence_df_based_on_fasta_library <- function(GI, Sequence, AA_in_protein, fixed_length, species){
  requireNamespace('stringr')
  requireNamespace('utils')
  # require(PhosMap)
  cat('Aligned sequence based on fasta library for motif enrichment anlysis.\n')
  # Read fasta file: gi_fasta
  PHOSPHATE_LIB_FASTA_FILE_PATH <- normalizePath(
    system.file(
      'extdata',
      'library', species, 'gi_fasta.txt',
      package = "PhosMap"
    ),
    mustWork = FALSE
  )
  if(!file.exists(PHOSPHATE_LIB_FASTA_FILE_PATH)){
    cat(PHOSPHATE_LIB_FASTA_FILE_PATH, ' -> ', 'No the file')
    stop('')
  }
  cat('Read fasta file of ', species, '.\n', sep = '')
  fasta_data <- utils::read.table(file=PHOSPHATE_LIB_FASTA_FILE_PATH, header=FALSE, sep="\t")
  border_limit <- floor(fixed_length/2)
  aligned_seq <- NULL
  GI_nrow <- length(GI)
  cat('Pre-align:', GI_nrow, 'phos-pepitdes.\n')
  cat('Fixed sequence length is ', fixed_length, '.\n', sep = '')
  cat('It needs few time.\n')
  for(i in seq_len(GI_nrow)){
    gi <- GI[i]
    aa_index <- AA_in_protein[i]
    loc_index <- as.numeric(stringr::str_split(aa_index, "[STY]", n = Inf, simplify = FALSE)[[1]])[2]
    index <- which(fasta_data[,1] == gi)
    if(length(index) > 0){
      refseq <- as.vector(fasta_data[index,2])
      refseq_len <- nchar(refseq)

      left_limit <- loc_index - border_limit
      right_limit <- loc_index + border_limit

      if(left_limit>=1 & right_limit>refseq_len){
        right_limit <- refseq_len
        truncated_seq <- stringr::str_sub(refseq, left_limit, right_limit)
        truncated_seq <- stringr::str_pad(truncated_seq, fixed_length, "right", pad = '_')
      }else if(left_limit<1 & right_limit<=refseq_len){
        left_limit <- 1
        truncated_seq <- stringr::str_sub(refseq, left_limit, right_limit)
        truncated_seq <- stringr::str_pad(truncated_seq, fixed_length, "left", pad = '_')
      }else if(left_limit<1 & right_limit>refseq_len){
        left_limit <- 1
        right_limit <- refseq_len
        truncated_seq <- stringr::str_sub(refseq, left_limit, right_limit)
        truncated_seq <- stringr::str_pad(truncated_seq, fixed_length, "both", pad = '_')
      }else{
        truncated_seq <- stringr::str_sub(refseq, left_limit, right_limit)
      }
    }else{
      truncated_seq <- NA
    }
    aligned_seq <- c(aligned_seq, truncated_seq)
    if(i %% 5000 == 0){
      cat('Aligned:', i, 'phos-pepitdes.\n')
    }
    if(i == GI_nrow){
      cat('Aligned:', i, 'phos-pepitdes.\n')
      cat('Finish OK! ^_^\n')
    }

  }
  cat('\n')
  aligned_sequence_df_based_on_fasta_library <- data.frame(GI, Sequence, AA_in_protein, aligned_seq)
  return(aligned_sequence_df_based_on_fasta_library)
}




