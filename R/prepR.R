#' umi4Cats prep
#'
#'@description
#'Prepare the fastq files for being analyse by the umi4Cats pipeline filtering reads with bait and
#'introducing the respective UMI in the file headers
#'
#'@usage
#'prepR(pathVenv, raw_dir, wk_dir, bait_seq, bait_pad, res_e)
#'
#'@inheritParams umi4CatsContacts
#'
#'\dontrun{
#'@examples
#'pathVenv = '/imppc/labs/lplab/share/marc/venv/umi4catsVenv
#'raw_dir = '/imppc/labs/lplab/share/marc/epimutations/processed/prove/raw'
#'wk_dir = '/imppc/labs/lplab/share/marc/epimutations/processed/prove/MLH1_ctrl_umi4cats_python'
#'threads = 1
#'bait_seq = 'ACCTAGAAGGATATGCGCTTGC'
#'bait_pad = 'GCGTTAGA'
#'res_e = "GATC"
#'
#'prepR(raw_dir, wk_dir, bait_seq, bait_pad, res_e)
#'}
#'
#'@exports

prepR <- function(pathVenv,
                  raw_dir,
                  wk_dir,
                  bait_seq,
                  bait_pad,
                  res_e){

  reticulate::use_virtualenv(pathVenv, required = T)
  py_functions <- system.file("python/umi4cats.py", package = "UMI4Cats")
  reticulate::source_python(py_functions)

  fastqmultx <- system.file("bin/fastq-multx/fastq-multx", package = "UMI4Cats")

  prep(raw_dir = raw_dir,
       wk_dir = wk_dir,
       bait_seq = bait_seq,
       bait_pad = bait_pad,
       res_e = res_e,
       fastqmultx = fastqmultx)

}


