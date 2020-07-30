dir_all <- downloadUMI4CexampleData(reduced = FALSE)
dir_red <- downloadUMI4CexampleData(reduced = TRUE)

test_that("Uncompressed example data folders exist", {
  expect_true(dir.exists(dir_all))
  expect_true(dir.exists(dir_red))

  expect_identical(dir(dir_all), c("CIITA", "ref_genome"))
  expect_identical(dir(dir_red), c("CIITA", "ref_genome"))
})

test_that("BFC is retrieved correctly", {
  expect_s4_class(.getCache(), "BiocFileCache")
})
