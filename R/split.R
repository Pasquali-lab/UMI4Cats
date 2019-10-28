#' Umi4c splitting
#'
#' This function splits the filtered reads by the restriction enzyme defined.
#' @inheritParams umi4CatsContacts
#' @export
#'
splitScript <- function(wk_dir,
                       cores,
                       re,
                       trimmomatic="trimmomatic"){
  # run prep script
  split_script <- system.file("exec/split.sh", package = "UMI4Cats")
  system(paste(split_script,
               "-w", wk_dir,
               "-c", cores,
               "-r", re,
               "-t", trimmomatic))

}
