#' Get viewpoint coordinates
#'
#' Finds the viewpoint coordinates for a given reference genome and sequence.
#'
#' @param bait_seq Character containing the bait primer sequence.
#' @param bait_pad Character containing the pad sequence (sequence between the
#' bait primer and the restriction enzyme sequence).
#' @param res_enz Character containing the restriction enzyme sequence.
#' @param ref_gen A BSgenome object of the reference genome.
#' @param sel_seqname A character with the chromosome name to focus the
#' search for the viewpoint sequence.
#' @return Creates a GRanges object containing the genomic position of the
#' viewpoint.
#' @examples
#' getViewpointCoordinates(
#'     bait_seq = "GGACAAGCTCCCTGCAACTCA",
#'     bait_pad = "GGACTTGCA",
#'     res_enz = "GATC",
#'     ref_gen = BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19,
#'     sel_seqname = "chr16" # Look only in chr16
#' )
#' @export
getViewpointCoordinates <- function(bait_seq,
    bait_pad,
    res_enz,
    ref_gen,
    sel_seqname=NULL) {
    viewpoint <- paste0(bait_seq, bait_pad, res_enz)

    # match viewpoint to genome
    if (is.null(sel_seqname)) {
        pos_viewpoint <- Biostrings::vmatchPattern(viewpoint,
                                                   ref_gen,
                                                   max.mismatch = 0
        )
    } else {
        ref_gen_sel <-  ref_gen[[sel_seqname]]
        pos_viewpoint <- Biostrings::matchPattern(viewpoint,
                                                  ref_gen_sel,
                                                  max.mismatch = 0
        )
        strand <- "+"
        
        if (length(pos_viewpoint) == 0) {
            pos_viewpoint <- Biostrings::matchPattern(Biostrings::reverseComplement(Biostrings::DNAString(viewpoint)),
                                                      ref_gen_sel,
                                                      max.mismatch = 0
            )
            strand <- "-"
        }
        pos_viewpoint <- GenomicRanges::GRanges(seqnames = sel_seqname,
                                                ranges = pos_viewpoint,
                                                strand = strand)
    }

    # only seqlevel of hit chromosome
    # seqlevels(pos_viewpoint) <- as.character(seqnames(pos_viewpoint))
    pos_viewpoint <- GenomeInfoDb::keepSeqlevels(
        pos_viewpoint,
        as.character(seqnames(pos_viewpoint))
    )

    return(pos_viewpoint)
}
