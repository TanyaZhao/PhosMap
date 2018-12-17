#' Get a data frame mapped GI number to Gene Symbol.
#'
#' This is an intermediate file and a dataframe with Gene Symbol is exported.
#' Based on a library file consisting of mapping relationships about Gene Symbol, GeneID and GI,
#' a new dataframe with Sequence, GI, Modification, Gene Symbol, Area and PSMs,is contructed.
#'
#' @param merge_df_with_phospho_peptides A dataframe consisting of IDs (Sequence_GI_Psite) and Area values.
#' @param species A string, the options is human and mouse, the default is human.
#'
#' @author Dongdong Zhan and Mengsha Tong
#'
#' @return A dataframe with Sequence, GI, Modification, Gene Symbol, Area values and PSMs
#' @export
#'
#' @examples
#' \dontrun{
#' combinated_df_with_mapped_gene_symbol = get_combined_data_frame_with_mapped_gene_symbol(
#'   merge_df_with_phospho_peptides
#' )
#' }

get_combined_data_frame_with_mapped_gene_symbol <- function(
  merge_df_with_phospho_peptides,
  species = 'human'
){
  # Read library file, map GI to Gene Symbol
  requireNamespace('utils')
  cat('\n The 5th step: write the data frame with symbols mapping to genes.')
  PHOSPHATE_LIB_MAPPING_FILE_PATH <- normalizePath(
    system.file(
      'extdata',
      'library', species, 'symbol_geneid_gi.txt',
      package = "PhosMap"
    ),
    mustWork = FALSE
  )
  if(!file.exists(PHOSPHATE_LIB_MAPPING_FILE_PATH)){
    cat(PHOSPHATE_LIB_MAPPING_FILE_PATH, ' -> ', 'No the file')
    stop('')
  }
  # Library file headers contains GeneSymbol, GeneID, GI
  symbol_geneid_gi <- utils::read.table(PHOSPHATE_LIB_MAPPING_FILE_PATH, header=T, sep="\t")

  cat('\n The 5th step is running.')
  # Split a string: sequenceID, accession, modification
  seq_gi_site_vector <- as.vector(merge_df_with_phospho_peptides$ID_of_seq_gi_site)
  Sequence <- apply(data.frame(seq_gi_site_vector), 1, function(x){
    strsplit(x, split="||", fixed = TRUE)[[1]][1]
  })
  GI <- apply(data.frame(seq_gi_site_vector), 1, function(x){
    strsplit(x, split="||", fixed = TRUE)[[1]][2]
  })
  Modification <- apply(data.frame(seq_gi_site_vector), 1, function(x){
    strsplit(x, split="||", fixed = TRUE)[[1]][3]
  })
  GeneSymbol <- apply(data.frame(GI), 1, function(x, MappingDf){
    gi.all <- strsplit(x, split=";", fixed = TRUE)[[1]]
    gi_mapping_symbol <- apply(data.frame(gi.all), 1, function(y, MappingDf){
      index_of_mapping <- which(MappingDf$GI==y)
      if(length(index_of_mapping)==0){
        return(NA)
      }else{
        symbol_of_mapping <- as.vector(MappingDf$GeneSymbol[index_of_mapping])
        return(symbol_of_mapping)
      }
    }, MappingDf = MappingDf)
    index_of_nonNA <- which(!is.na(gi_mapping_symbol))
    if(length(index_of_nonNA)>0){
      gi_mapping_symbol <- gi_mapping_symbol[index_of_nonNA]
      gi_mapping_symbol <- unique(gi_mapping_symbol)
      gi_mapping_symbol <- paste(gi_mapping_symbol, collapse=";")
      return(gi_mapping_symbol)
    }else{
      return(NA)
    }
  }, MappingDf = symbol_geneid_gi)

  # sequenceID, accession, symbol, modification, quantification_value_in_experiment
  df_of_combination <- data.frame(Sequence, GI, Modification, GeneSymbol, merge_df_with_phospho_peptides[,-1]) # delete first column
  index_of_NonNA <- which(!is.na(GeneSymbol))
  df_of_combination <- df_of_combination[index_of_NonNA,]
  cat('\n The 5th step is over ^_^.')
  cat('\n The 5th step: write the data frame with symbols mapping to genes.')
  return(df_of_combination)
}
