path <- digestGenome(
  res_enz = "GATC",
  cut_pos = 0,
  name_RE = "dpnII",
  sel_chr = "chrY",
  ref_gen = BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19,
  out_path = file.path(tempdir(), "digested_genome/")
)

test_that("digestGenome directory and files are created", {
  expect_true(dir.exists(path))
  expect_true(dir(path)>0)
  expect_true(all(grepl(".rda", list.files(path))))
})

test_that("digestGenome output is correct", {
  load(list.files(path, full.names=TRUE))

  expect_equal(ls(), "digested_genome_gr")
  expect_equal(length(digested_genome_gr), 60005)

  seqs <- getSeq(BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19, digested_genome_gr[-c(1,length(digested_genome_gr))])
  expect_true(all(grepl("^GATC", seqs)))
})
