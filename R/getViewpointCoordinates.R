#' Get viewpoint coordinates
#'
#' Finds the viewpoint coordinates for a given reference genome.
#'
#' @param bait_seq Character containing the bait primer sequence.
#' @param bait_pad Character containing the pad sequence (sequence between the bait primer and the restriction enzyme sequence).
#' @param res_enz Character containing the restriction enzyme sequence.
#' @param ref_gen A BSgenome object of the reference genome.
#' @return Creates a GRanges object containing the genomic position of the viewpoint.
#' @examples
#' \dontrun{
#' getViewpointCoordinates(bait_seq="CCCAAATCGCCCAGACCAG",
#'               bait_pad="GCGCG",
#'               res_enz="GATC",
#'               ref_gen=BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19)
#' }
#' @export
getViewpointCoordinates <- function(bait_seq,
                                    bait_pad,
                                    res_enz,
                                    ref_gen) {

  viewpoint <-  paste0(bait_seq, bait_pad, res_enz)

  # match viewpoint to genome
  pos_viewpoint <- Biostrings::vmatchPattern(viewpoint, ref_gen, max.mismatch = 0)

  # only seqlevel of hit chromosome
  seqlevels(pos_viewpoint) <- as.character(seqnames(pos_viewpoint))

  return(pos_viewpoint)
}


