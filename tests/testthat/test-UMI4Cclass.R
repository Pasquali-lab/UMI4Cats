umi <- makeUMI4Cexample()

enh <- GRanges(c(
  "chr16:10925006-10928900",
  "chr16:11102721-11103700"
))

# Perform differential test
umi_dif <- fisherUMI4C(umi,
                      query_regions = enh,
                        filter_low = 20,
                        resize = 5e3
)

test_that("Example UMI4C object is created successfully", {
  expect_s4_class(umi, "UMI4C")
  expect_equal(length(umi), 2838)
  expect_equal(colnames(umi), c("ctrl", "cyt"))
})

test_that("UMI4C methods work correctly", {

  # Bait method
  expect_s4_class(bait(umi), "GRanges")

  # dgram method
  expect_s4_class(dgram(umi), "SimpleList")
  expect_true(length(dgram(umi)$ctrl)==length(dgram(umi)$cyt))

  # trend method
  expect_equal(class(trend(umi)), "data.frame")

  # results method
  expect_error(results(umi))
  expect_s4_class(results(umi_dif), "GRanges")
})
