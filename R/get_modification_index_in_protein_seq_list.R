#' Get index of modifications in protein sequence.
#' @param id_data_only_peptide2gi a data frame for peptides with protein gi.
#' @param gi_fasta a fasta data for a specific species.
#'
#' @author Dongdong Zhan and Mengsha Tong
#'
#' @return A vector for index of modifications in protein sequence.
#' @export
#'
get_modification_index_in_protein_seq_list <- function(id_data_only_peptide2gi, gi_fasta){
  # 1
  # Get modification index in protein sequence.
  cat('\n', 'Get modification index in protein sequence.')
  id_data_only_peptide2gi_row <- nrow(id_data_only_peptide2gi)
  modification_index_in_protein_seq_list <- list()
  for(i in seq_len(id_data_only_peptide2gi_row)){
    peptide_seq <- as.vector(id_data_only_peptide2gi$Sequence[i])
    peptide_gi <- as.vector(id_data_only_peptide2gi$GI[i])
    modification_index_in_peptide_seq <- unlist(gregexpr("[a-z]", peptide_seq))
    protein_seq <- as.vector(gi_fasta$Sequence[which(gi_fasta$GI==peptide_gi)])
    first_index_of_peptide2protein <- unlist(gregexpr(toupper(peptide_seq), protein_seq))
    modification_index_in_protein_seq <- NULL
    for(elemt in first_index_of_peptide2protein){
      tmp_modification_index_in_protein_seq <- elemt + modification_index_in_peptide_seq -1
      modification_index_in_protein_seq <- c(modification_index_in_protein_seq,
                                             tmp_modification_index_in_protein_seq)
    }
    modification_index_in_protein_seq_list[[i]] <- modification_index_in_protein_seq
    if(i%%500==0 | i==id_data_only_peptide2gi_row ){
      cat('\n completed: ', i, '/', id_data_only_peptide2gi_row)
    }
  }
  return(modification_index_in_protein_seq_list)
}
