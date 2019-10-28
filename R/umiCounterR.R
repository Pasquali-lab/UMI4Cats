#' umi4Cats umiCounter
#' 
#'@description
#'Count number of UMIs for every ligation taking ligation position and UMI mismatches into acount. It generates a tsv file 
#'for creating the 4c profile. 
#'  
#'@usage
#'umiCounterR(py_umi4cats, raw_dir, wk_dir, bait_seq, bait_pad, res_e, fastqmultx)
#'  
#'@inheritParams umi4CatsContacts
#'  
#'@example 
#'py_umi4cats = '/imppc/labs/lplab/share/marc/repos/umi4cats/python/umi4catsBuilder.py'
#'wk_dir = '/imppc/labs/lplab/share/marc/epimutations/processed/prove/MLH1_ctrl_umi4cats_python'
#'bait_seq = 'ACCTAGAAGGATATGCGCTTGC'
#'bait_pad = 'GCGTTAGA'
#'res_e = "GATC"
#'bowtie2 = 'bowtie2'
#'ref_gen = '/imppc/labs/lplab/share/marc/refgen/hg19/hg19.fa'
#'samtools = 'samtools'
#'genomic_track = '/imppc/labs/lplab/share/marc/epimutations/processed/genomicTracks/genomic_tracks_hg19/dpnII_genomicTrack'
#'  
#'umiCounterR(py_umi4cats, wk_dir, bowtie2, ref_gen, samtools, genomicTrack, bait_seq, bait_pad,res_e)
#'
#'@export

umiCounterR <- function(py_umi4cats,
                        wk_dir,
                        bowtie2,
                        ref_gen,
                        samtools,
                        genomic_track,
                        bait_seq,
                        bait_pad,
                        res_e){
  
  library(reticulate)
  use_python(py_umi4cats, required = T)
  py_functions <- system.file("python/umi4cats.py", package = "UMI4Cats")
  source_python(py_functions)
  
  umiCounter(wk_dir,
             bowtie2,
             ref_gen,
             samtools,
             genomic_track,
             bait_seq,
             bait_pad,
             res_e)
  
}


