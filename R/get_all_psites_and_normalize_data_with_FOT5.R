#' To normalize data and filter data only including phosphorylation sites.
#'
#' @param data_frame A data frame containing IDs and quantification values merged from multi-experiments as input.
#' @param experiment_code_file_path A file path of storing experiment codes as input. The experiment codes are required to keep pace with column names of Value.
#'
#' @return A data frame after normalization and filtering (x 1e5)
#' @export
#'
#' @examples
#' \dontrun{
#' ptypes_df = get_normalized_data_FOT5(data_frame, experiment_code_file_path)
#' }

get_all_psites_and_normalize_data_with_FOT5 <- function(data_frame, experiment_code_file_path){
  requireNamespace('utils')
  cat('\n The 7th step: Normalize data and filter data only including phosphorylation site.')
  experiment_code <- utils::read.table(experiment_code_file_path, header = TRUE, sep = '\t', stringsAsFactors = NA)
  experiment_code <- as.vector(unlist(experiment_code$Experiment_Code))
  data_frame_colnames <- colnames(data_frame)

  cat('\n The 7th step is running.')
  summary_df_ID_Info <- data_frame[, seq_len(6)]
  summary_df_Value <- data_frame[, -(seq_len(6))]
  Value_Area <- summary_df_Value
  Value_FOT5 <- Value_Area
  Value_FOT5_col <- ncol(Value_FOT5)
  for(i in seq_len(Value_FOT5_col)){
    x <- Value_Area[,i]
    valid_index <- which(x>0)
    valid_x <- x[valid_index]
    valid_x_sum <- sum(valid_x)
    valid_x_FOT5 <- valid_x/valid_x_sum*1e5
    Value_FOT5[valid_index,i] <- valid_x_FOT5
  }
  summary_df_Value_Fot <- as.matrix( Value_FOT5)

  cat('\n Imputation with the next order of magnitude of the minimum except for zero.')
  index_of_zero <- which(summary_df_Value_Fot==0)
  min_value_of_non_zero <- min(summary_df_Value_Fot[-index_of_zero])
  summary_df_Value_Fot[index_of_zero] <- min_value_of_non_zero*0.1
  summary_df_Value_Fot_With_ID <- data.frame(summary_df_ID_Info, summary_df_Value_Fot)
  summary_df_Value_Fot_With_ID$AA_in_protein <- toupper(summary_df_Value_Fot_With_ID$AA_in_protein)

  cat('\n Filtering data only including S/T/Y modifications.')
  ptypes <- c('S', 'T', 'Y')
  index_of_AA_in_protein <- apply(data.frame(summary_df_Value_Fot_With_ID$AA_in_protein), 1, function(x){
    if(grepl('S', x) | grepl('T', x) | grepl('Y', x)){
      return(TRUE)
    }else{
      return(FALSE)
    }
  })
  index_of_ptypes <- which(index_of_AA_in_protein)
  if(length(index_of_ptypes)>0){
    ptypes_df <- summary_df_Value_Fot_With_ID[index_of_ptypes,]
  }else{
    stop('No phosphorylation data')
  }
  cat('\n The 7th step is over ^_^.')
  return(ptypes_df)
}
