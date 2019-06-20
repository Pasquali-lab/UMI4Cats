filterRaw <- function(input,
                      bait_exclusion=3e3,
                      bait_coordinates,
                      bait_upstream=5e5,
                      bait_downstream=5e5) {
  ## Order input by coordinates
  input <- input[order(input[,1]),]

  ## Remove region around bait
  bait_exp <- GenomicRanges::resize(bait_coordinates, fix="center", width=bait_exclusion)
  l <- input[,1] <= GenomicRanges::start(bait_exp) |
    input[,1] >= GenomicRanges::end(bait_exp)
  input <- input[l,]

  ## Remove regions outside region
  region <- regioneR::extendRegions(bait_coordinates,
                                    extend.start=bait_upstream,
                                    extend.end=bait_downstream)

  l <- input[,1] >= GenomicRanges::start(region) &
    input[,1] <= GenomicRanges::end(region)
  input <- input[l,]

  ## Create raw output
  raw <- as.data.frame(input[,1:2])
  colnames(raw) <- c("coord", "UMIs")

  return(raw)
}
