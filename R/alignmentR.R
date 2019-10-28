#' umi4Cats alignment
#' 
#'@description
#' Align splitted UMI reads agains a pre-defined reference genome using Bowtie2
#'  
#'@usage
#' alignmentR(py_umi4cats, raw_dir, wk_dir, bait_seq, bait_pad, res_e, fastqmultx)
#'  
#'@inheritParams umi4CatsContacts
#'  
#'@example
#'py_umi4cats = '/imppc/labs/lplab/share/marc/repos/umi4cats/python/umi4catsBuilder.py'
#'wk_dir = '/imppc/labs/lplab/share/marc/epimutations/processed/prove/MLH1_ctrl_umi4cats_python'
#'threads = 1
#'bait_seq = 'ACCTAGAAGGATATGCGCTTGC'
#'bait_pad = 'GCGTTAGA'
#'res_e = "GATC"
#'cut_pos = 0
#'bowtie2 = 'bowtie2'
#'ref_gen = '/imppc/labs/lplab/share/marc/refgen/hg19/hg19.fa'
#'samtools = 'samtools'
#'  
#'alignmentR(py_umi4cats, wk_dir, threads, bowtie2, ref_gen, samtools, bait_seq, bait_pad, res_e)
#'  
#'@export

alignmentR <- function(py_umi4cats,
                       wk_dir,
                        threads,
                        bowtie2,
                        ref_gen,
                        samtools,
                        bait_seq,
                        bait_pad,
                        res_e){
  
  library(reticulate)
  use_python(py_umi4cats, required = T)
  py_functions <- system.file("python/umi4cats.py", package = "UMI4Cats")
  source_python(py_functions)
  
  alignment(wk_dir,
            threads,
            bowtie2,
            ref_gen,
            samtools,
            bait_seq,
            bait_pad,
            res_e)
  
}
