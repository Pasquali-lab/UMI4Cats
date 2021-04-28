coord <- getViewpointCoordinates(
      bait_seq = "GCGTGTAGCCCCCAGTCACGTTCCGA",
      bait_pad = "TGCTTGCTC",
      res_enz = "GATC",
      ref_gen = BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19,
      sel_seqname = "chrY" # Look only in chr16
  )

test_that("Output coordinates are correct", {
  expect_s4_class(coord, "GRanges")
  expect_equal(as.character(coord), "chrY:12560-12598:+")
  expect_equal(as.character(seqnames(coord)), "chrY")
})
