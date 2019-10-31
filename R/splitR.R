#' umi4Cats splitting
#'
#'@description
#' Split the prepared reads using the pre-defined restrition enzyme information
#'
#'@usage
#'splitR(pathVenv, wk_dir, res_e, cut_pos)
#'
#'@inheritParams umi4CatsContacts
#'
#'\dontrun{
#'@examples
#'pathVenv = '/imppc/labs/lplab/share/marc/venv/umi4catsVenv
#'wk_dir = '/imppc/labs/lplab/share/marc/epimutations/processed/prove/MLH1_ctrl_umi4cats_python'
#'res_e = "GATC"
#'cut_pos = 0
#'
#'splitR(wk_dir, res_e, cut_pos)
#'}
#'
#'@export

splitR <- function(pathVenv,
                   wk_dir,
                  res_e,
                  cut_pos){


  reticulate::use_virtualenv(pathVenv, required = T)
  py_functions <- system.file("python/umi4cats.py", package = "UMI4Cats")
  reticulate::source_python(py_functions)

  cut_pos <- as.character(cut_pos)

  split(wk_dir = wk_dir,
        res_e = res_e,
        cut_pos = cut_pos)

}
