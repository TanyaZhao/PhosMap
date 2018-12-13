#' Generate peptide identification files with psites scores.
#'
#' Based on mascot txt files with psites and peptide identification files downloaded from Firmiana, the file including phosphorylation modifications is generated.
#'
#'
#' @param mascot_txt_dir A folder containing identification xml files with psites scores as input.
#' @param firmiana_peptide_dir A folder containing identification txt file downloaded from Firmiana as input.
#' @param psites_score_dir A folder used for saving files of peptide identification file with psites scores
#'
#' @author Dongdong Zhan and Mengsha Tong
#'
#' @return A series of output files saved in the psites_score_dir
#'
#' @examples
#' \dontrun{
#' generate_psites_score_file(mascot_txt_dir, firmiana_peptide_dir, psites_score_dir)
#' }

generate_psites_score_file <- function(mascot_txt_dir, firmiana_peptide_dir, psites_score_dir){
  mascot_txt_dir_paths <- list.dirs(mascot_txt_dir)
  mascot_txt_dir_paths <- mascot_txt_dir_paths[-1]
  mascot_txt_dir_paths.expNames <- apply(data.frame(mascot_txt_dir_paths), 1, function(x){
    x <- strsplit(x, split = '/')[[1]]
    x[length(x)]
  })

  firmiana_txt_file_names <- list.files(firmiana_peptide_dir)
  firmiana_txt_file_names.expNames <- apply(data.frame(firmiana_txt_file_names), 1, function(x){
    x <- strsplit(x, split = '_')[[1]]
    x[1]
  })
  mascot_txt_dir_paths.len <- length(mascot_txt_dir_paths)
  cat('\n Total file: ', mascot_txt_dir_paths.len)
  cat('\n It will take a little while...')
  for(i in seq_len(mascot_txt_dir_paths.len)){
    mascot_txt_dir_path <- mascot_txt_dir_paths[i]
    mascot_txt_dir_path.expName <- mascot_txt_dir_paths.expNames[i]
    mascot_txt_dir_path_expName_path <- list.files(mascot_txt_dir_path)
    mascot_txt_dir_path_expName_path <- normalizePath(paste(mascot_txt_dir_path, mascot_txt_dir_path_expName_path, sep = '/'))
    match.index <- match(mascot_txt_dir_path.expName, firmiana_txt_file_names.expNames)
    firmiana_txt_file_name <- firmiana_txt_file_names[match.index]
    firmiana_peptide_dir_path_expName_path <- normalizePath(paste(firmiana_peptide_dir, firmiana_txt_file_name, sep = '/'))

    outputName <- normalizePath(paste(psites_score_dir, '/', mascot_txt_dir_path.expName, '_psites_score.csv', sep = ''), mustWork = F)
    expName <- mascot_txt_dir_path.expName
    write_csv_pep_seq_conf(expName, outputName, mascot_txt_dir_path_expName_path, firmiana_peptide_dir_path_expName_path)

    cat('\n Completed file: ', i, '/', mascot_txt_dir_paths.len)

  }

}





