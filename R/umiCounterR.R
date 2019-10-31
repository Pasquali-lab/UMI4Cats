#' umi4Cats umiCounter
#'
#'@description
#'Count number of UMIs for every ligation taking ligation position and UMI mismatches into acount. It generates a tsv file
#'for creating the 4c profile.
#'
#'@usage
#'umiCounterR(pathVenv, wk_dir, ref_gen, genomic_track, bait_seq, bait_pad, res_e)

#'
#'@inheritParams umi4CatsContacts
#'
#'@examples
#'\dontrun{
#'pathVenv = '/imppc/labs/lplab/share/marc/venv/umi4catsVenv'
#'wk_dir = '/imppc/labs/lplab/share/marc/epimutations/processed/prove/MLH1_ctrl_umi4cats_python'
#'bait_seq = 'ACCTAGAAGGATATGCGCTTGC'
#'bait_pad = 'GCGTTAGA'
#'res_e = "GATC"
#'ref_gen = '/imppc/labs/lplab/share/marc/refgen/hg19/hg19.fa'
#'genomic_track = '/imppc/labs/lplab/share/marc/epimutations/processed/genomicTracks/genomic_tracks_hg19/dpnII_genomicTrack'
#'
# umiCounterR(pathVenv = pathVenv,
#             wk_dir = wk_dir,
#             ref_gen = ref_gen,
#             genomic_track = genomic_track,
#             bait_seq = bait_seq,
#             bait_pad = bait_pad,
#             res_e = res_e)
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

}


