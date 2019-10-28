#' umi4Cats splitting
#' 
#'@description
#' Split the prepared reads using the pre-defined restrition enzyme information
#'  
#'@usage 
#'splitR(py_umi4cats, wk_dir, res_e, cut_pos)
#'  
#'@inheritParams umi4CatsContacts
#'  
#'@example
#' py_umi4cats = '/imppc/labs/lplab/share/marc/repos/umi4cats/python/umi4catsBuilder.py'
#'wk_dir = '/imppc/labs/lplab/share/marc/epimutations/processed/prove/MLH1_ctrl_umi4cats_python'
#'res_e = "GATC"
#'cut_pos = 0 
#'splitR(py_umi4cats, wk_dir, res_e, cut_pos)
#'@export

splitR <- function(py_umi4cats,
                   wk_dir,
                  res_e,
                  cut_pos){
  

  library(reticulate)
  use_python(py_umi4cats, required = T)
  py_functions <- system.file("python/umi4cats.py", package = "UMI4Cats")
  source_python(py_functions)
  
  cut_pos <- as.character(cut_pos)
  
  split(wk_dir,
        res_e,
        cut_pos)
  
}
