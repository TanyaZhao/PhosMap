
#' Assign psites to protein sequence.
#'
#' Construct the data frame with unique phosphorylation site for each protein sequence and eliminate redundancy.
#'
#'
#' @param combinated_df_with_mapped_gene_symbol A dataframe with Sequence, GI, Modification, Gene Symbol, Area and PSMs as input.
#' @param species A string, the options is human and mouse, the default is human.
#'
#' @author Dongdong Zhan and Mengsha Tong
#'
#' @return A dataframe that all redundant psites are assigned to protein sequence.
#' @export
#'
#' @examples
#' \dontrun{
#' summary_df_of_unique_proteins_with_sites = get_summary_df_of_unique_proteins_with_sites(
#'   combinated_df_with_mapped_gene_symbol
#' )
#' }
get_summary_df_of_unique_proteins_with_sites <- function(
  combinated_df_with_mapped_gene_symbol,
  species = 'human'
){
  requireNamespace('utils')
  # unique phosphorylation sites
  cat('\n The 6th step: construct the data frame with unique phosphorylation site for each protein sequence.')
  # Read fasta file: gi_fasta
  PHOSPHATE_LIB_FASTA_FILE_PATH <- normalizePath(
    system.file(
      'extdata',
      'library', species, 'gi_fasta.txt',
      package = "PhosMap"
    ),
    mustWork = F
  )
  if(!file.exists(PHOSPHATE_LIB_FASTA_FILE_PATH)){
    cat(PHOSPHATE_LIB_FASTA_FILE_PATH, ' -> ', 'No the file')
    stop('')
  }
  gi_fasta <- utils::read.table(file=PHOSPHATE_LIB_FASTA_FILE_PATH, header=F, sep="\t")
  colnames(gi_fasta) <- c('GI', 'Sequence')


  id_data <- combinated_df_with_mapped_gene_symbol

  # Keep peptides assigned to unique protein
  id_data_only_peptide2gi <- id_data[which(!grepl(';', as.vector(id_data$GI))),]

  # Determine locations of the psites each peptide mapped to protein squence.
  modification_index_in_protein_seq_list <- get_modification_index_in_protein_seq_list(id_data_only_peptide2gi,
                                                                                       gi_fasta)

  proteins_in_id_data_only_peptide2gi <- as.vector(id_data_only_peptide2gi$GI)
  sequences_in_id_data_only_peptide2gi <- as.vector(id_data_only_peptide2gi$Sequence)
  value_in_id_data_only_peptide2gi <- id_data_only_peptide2gi[, -c(1:4)]

  unique_proteins <- unique(proteins_in_id_data_only_peptide2gi)
  unique_protein_count <- length(unique_proteins)

  # Show psites and modifications of one protein, merge the values with the same modification type.
  cat('\n', 'Map phosphorylation sites to protein sequence and eliminate redundancy.')
  system.time({
    summary_df_of_unique_proteins_with_sites <- c()
    for(i in seq_len(unique_protein_count)){

      df_with_AAs_i <- get_df_with_AAs_i(unique_proteins,
                                        i,
                                        id_data_only_peptide2gi,
                                        proteins_in_id_data_only_peptide2gi,
                                        sequences_in_id_data_only_peptide2gi,
                                        modification_index_in_protein_seq_list)


      summary_df_of_unique_protein_with_sites <- get_df_with_AAs_i_of_reducing_redundancy(df_with_AAs_i)


      summary_df_of_unique_proteins_with_sites <- rbind(summary_df_of_unique_proteins_with_sites, summary_df_of_unique_protein_with_sites)

      if(i%%500==0 | i == unique_protein_count){
        cat('\n completed: ', i, '/', unique_protein_count)
      }
    }
  })
  summary_df_of_unique_proteins_with_sites_rownames <- paste(as.vector(summary_df_of_unique_proteins_with_sites$GI),
                                                             as.vector(summary_df_of_unique_proteins_with_sites$AA_in_protein),
                                                             sep = '_')
  rownames(summary_df_of_unique_proteins_with_sites) <- summary_df_of_unique_proteins_with_sites_rownames
  summary_df_of_unique_proteins_with_sites_colnames <- colnames(summary_df_of_unique_proteins_with_sites)
  index_of_PSMs <- which(grepl('_PSMs', summary_df_of_unique_proteins_with_sites_colnames))
  if(length(index_of_PSMs)>0){
    summary_df_of_unique_proteins_with_sites <- summary_df_of_unique_proteins_with_sites[,-index_of_PSMs]
  }
  summary_df_of_unique_proteins_with_sites$GeneSymbol <- apply(data.frame(summary_df_of_unique_proteins_with_sites$GeneSymbol),
                                                               1,
                                                               function(x){
    if(grepl('||', x)){
      x <- as.vector(x)
      x <- strsplit(x, split = '||', fixed = T)
      x[[1]][1]
    }
  })
  cat('\n The 6th step: construct over.')
  return(summary_df_of_unique_proteins_with_sites)
}






