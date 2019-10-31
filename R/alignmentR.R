#' umi4Cats alignment
#'
#'@description
#' Align splitted UMI reads agains a pre-defined reference genome using Bowtie2
#'
#'@usage
#'alignmentR(pathVenv, wk_dir, threads, ref_gen, bait_seq, bait_pad, res_e)
#'
#'@inheritParams umi4CatsContacts
#'
#'@examples
#'\dontrun{
#'pathVenv = '/imppc/labs/lplab/share/marc/venv/umi4catsVenv'
#'wk_dir = '/imppc/labs/lplab/share/marc/epimutations/processed/prove/MLH1_ctrl_umi4cats_python'
#'threads = 1
#'bait_seq = 'ACCTAGAAGGATATGCGCTTGC'
#'bait_pad = 'GCGTTAGA'
#'res_e = "GATC"
#'cut_pos = 0
#'ref_gen = '/imppc/labs/lplab/share/marc/refgen/hg19/hg19.fa'
#'
#'alignmentR(wk_dir, threads, ref_gen, bait_seq, bait_pad, res_e)
#'}
#'
#'@export

alignmentR <- function(pathVenv,
                       wk_dir,
                        threads,
                        ref_gen,
                        bait_seq,
                        bait_pad,
                        res_e){

  reticulate::use_virtualenv(pathVenv, required = T)
  py_functions <- system.file("python/umi4cats.py", package = "UMI4Cats")
  reticulate::source_python(py_functions)

  bowtie2 <- 'bowtie2'
  samtools <- 'samtools'

  alignment(wk_dir = wk_dir,
            threads = threads,
            bowtie2 = bowtie2,
            ref_gen = ref_gen,
            samtools = samtools,
            bait_seq = bait_seq,
            bait_pad = bait_pad,
            res_e = res_e)

}
