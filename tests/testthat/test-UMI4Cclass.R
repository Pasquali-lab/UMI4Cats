data("ex_ciita_umi4c")
ex_ciita_umi4c <- addGrouping(ex_ciita_umi4c, grouping="condition")
umi <- ex_ciita_umi4c

enh <- GRanges(c(
  "chr16:10925006-10928900",
  "chr16:11102721-11103700"
))

# Perform differential test
umi_dif <- fisherUMI4C(umi,
                      grouping = "condition",
                      query_regions = enh,
                      filter_low = 20,
                      resize = 5e3
)

test_that("Example UMI4C object is created successfully", {
  expect_s4_class(umi, "UMI4C")
  expect_equal(colnames(groupsUMI4C(umi)[["condition"]]), c("ctrl", "cyt"))
})

test_that("UMI4C methods work correctly", {

  # Bait method
  expect_s4_class(bait(umi), "GRanges")

  # dgram method
  expect_s4_class(dgram(umi), "SimpleList")
  expect_true(length(dgram(umi)[[1]])==length(dgram(umi)[[2]]))

  # trend method
  expect_equal(class(trend(umi)), "data.frame")

  # results method
  expect_error(resultsUMI4C(umi))
  expect_s4_class(resultsUMI4C(umi_dif), "GRanges")
})
