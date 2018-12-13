#' Write data to specific direction with CSV format.
#'
#'
#' @param expName a string for experiment name as input.
#' @param outputName a string for experiment name as output.
#' @param mascotfileNames a vector for storing mascot file names.
#' @param refFileName a string for reference file name.
#'
#' @author Dongdong Zhan and Mengsha Tong
#'
#' @export


write_csv_pep_seq_conf <- function(expName, outputName, mascotfileNames, refFileName){
  requireNamespace('utils')
  fileNames.len <- length(mascotfileNames)
  fileData.list <- list()
  for(i in 1:fileNames.len){
    mascotfileName <- mascotfileNames[i]
    df <- get_filtered_df(mascotfileName, refFileName)
    fileData.list[[i]] <- df
  }

  mergeMat <- subset(fileData.list[[1]], select = c('pep_seq', 'pep_var_mod_conf'))
  if(fileNames.len>1){
    for(i in 2:fileNames.len){
      tmp.mergeMat <- subset(fileData.list[[i]], select = c('pep_seq', 'pep_var_mod_conf'))
      mergeMat <- merge(mergeMat, tmp.mergeMat, by='pep_seq', all = T)
    }
  }

  mergeMat_ionscore <- subset(fileData.list[[1]], select = c('pep_seq', 'pep_score'))
  if(fileNames.len>1){
    for(i in 2:fileNames.len){
      tmp.mergeMat_ionscore <- subset(fileData.list[[i]], select = c('pep_seq', 'pep_score'))
      mergeMat_ionscore <- merge(mergeMat_ionscore, tmp.mergeMat_ionscore, by='pep_seq', all = T)
    }
  }

  mergeMat.row <- nrow(mergeMat)
  pep_seq <- as.vector(mergeMat$pep_seq)
  pep_var_mod_conf <- NULL
  for(i in seq_len(mergeMat.row)){
    score <- as.vector(unlist(mergeMat_ionscore[i,-1]))
    conf <- as.vector(unlist(mergeMat[i,-1]))
    index <- which.max(score)
    conf <- conf[index]
    pep_var_mod_conf <- c(pep_var_mod_conf, conf)
  }

  df1 <- data.frame(pep_seq, pep_var_mod_conf)
  utils::write.csv(df1, outputName, row.names = F)

}



