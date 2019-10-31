#' umi4Cats umiCounter
#'
#'@description
#'Count number of UMIs for every ligation taking ligation position and UMI mismatches into acount. It generates a tsv file
#'for creating the 4c profile.
#'
#'@usage
#'umiCounterR(pathVenv, raw_dir, wk_dir, bait_seq, bait_pad, res_e, fastqmultx)
#'
#'@inheritParams umi4CatsContacts
#'
#'\dontrun{
#'@examples
#'pathVenv = '/imppc/labs/lplab/share/marc/venv/umi4catsVenv
#'wk_dir = '/imppc/labs/lplab/share/marc/epimutations/processed/prove/MLH1_ctrl_umi4cats_python'
#'bait_seq = 'ACCTAGAAGGATATGCGCTTGC'
#'bait_pad = 'GCGTTAGA'
#'res_e = "GATC"
#'ref_gen = '/imppc/labs/lplab/share/marc/refgen/hg19/hg19.fa'
#'genomic_track = '/imppc/labs/lplab/share/marc/epimutations/processed/genomicTracks/genomic_tracks_hg19/dpnII_genomicTrack'
#'
#'umiCounterR(wk_dir, ref_gen, genomicTrack, bait_seq, bait_pad,res_e)
#'
#'}
#'
#'@export

umiCounterR <- function(pathVenv,
                        wk_dir,
                        ref_gen,
                        genomic_track,
                        bait_seq,
                        bait_pad,
                        res_e){

  reticulate::use_virtualenv(pathVenv, required = T)
  py_functions <- system.file("python/umi4cats.py", package = "UMI4Cats")
  reticulate::source_python(py_functions)

  bowtie2 <- 'bowtie2'
  samtools <- 'samtools'

  umiCounter(wk_dir = wk_dir,
             bowtie2 = bowtie2,
             ref_gen = ref_gen,
             samtools = samtools,
             genomic_track = genomic_track,
             bait_seq = bait_seq,
             bait_pad = bait_pad,
             res_e = res_e)


  # umi4counts <- list.files(file.path(wk_dir, 'rst'),
  #                          pattern = "_umi_counts.tsv",
  #                          full.names = T)
  #
  #
  # umi4counts <- read.table(umi4count)

}


