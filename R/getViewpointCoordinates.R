getViewpointCoordinates <- function(bait_seq,
                                    bait_pad,
                                    res_enz,
                                    ref_gen) {
  viewpoint <-  paste0(bait_seq, bait_pad, res_enz)

  # TODO: bowtie index autogeneration if not exist? set automatic bowtie or define path?
  bowtie_index <- gsub('\\.fa$', '', ref_gen)

  # TODO: Define viewpoint in main function to avoid doing several searches in each function.
  view_point_pos <- system(paste(system.file(package="Rbowtie2", "bowtie2-align-s"),
                                 '--quiet',
                                 '-N 0',
                                 '-x', bowtie_index,
                                 '-c', viewpoint),
                           intern = T)


  view_point_pos <- tail(view_point_pos, n = 1)
  view_point_pos <- unlist(strsplit(view_point_pos, "\t"))
  pos_chr <- view_point_pos[3]
  pos_start <- as.numeric(view_point_pos[4])
  pos_end <- pos_start + nchar(viewpoint) - nchar(res_enz)
  pos_viewpoint <- GenomicRanges::GRanges(seqnames = pos_chr,
                                          ranges = IRanges::IRanges(start=pos_start,
                                                                   end=pos_end)
                                          )

  return(pos_viewpoint)
}
