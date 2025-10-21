window <- GRanges("chr16:11348649-11349648")
gene_anno <- createGeneAnnotation(
    window = window,
    TxDb = TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene,
    OrgDb = "org.Hs.eg.db"
)

test_that("Output of createGeneAnnotation is correct", {
  expect_s4_class(gene_anno, "GRanges")
  expect_equal(unique(strand(gene_anno)), factor(c("+", "-"), levels=c("+", "-", "*")))
  expect_equal(unique(gene_anno$type), c("GENE", "EXON"))
})

