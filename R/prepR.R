#' umi4Cats prep
#' 
#'@description 
#'Prepare the fastq files for being analyse by the umi4Cats pipeline filtering reads with bait and
#'introducing the respective UMI in the file headers 
#'  
#'@usage 
#'prepR(py_umi4cats, raw_dir, wk_dir, bait_seq, bait_pad, res_e, fastqmultx)
#'  
#'@inheritParams umi4CatsContacts
#'  
#'@examples
#'py_umi4cats = '/imppc/labs/lplab/share/marc/repos/umi4cats/python/umi4catsBuilder.py'
#'raw_dir = '/imppc/labs/lplab/share/marc/epimutations/processed/prove/raw'
#'wk_dir = '/imppc/labs/lplab/share/marc/epimutations/processed/prove/MLH1_ctrl_umi4cats_python'
#'threads = 1
#'bait_seq = 'ACCTAGAAGGATATGCGCTTGC'
#'bait_pad = 'GCGTTAGA'
#'res_e = "GATC"
#'fastqmultx = 'fastq-multx'
#'  
#'prepR(py_umi4cats, raw_dir, wk_dir, bait_seq, bait_pad, res_e, fastqmultx)
#'  
#'@export

prepR <- function(py_umi4cats,
                  raw_dir,
                  wk_dir,
                  bait_seq,
                  bait_pad,
                  res_e,
                  fastqmultx){
  
  library(reticulate)
  use_python(py_umi4cats, required = T)
  py_functions <- system.file("python/umi4cats.py", package = "UMI4Cats")
  source_python(py_functions)
  
  
  prep(raw_dir,
       wk_dir,
       bait_seq,
       bait_pad,
       res_e,
       fastqmultx)
  
}


  