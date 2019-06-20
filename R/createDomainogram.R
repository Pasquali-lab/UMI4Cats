#' Create Domainogram
#'
#' Using as input the raw UMIs, this package creates a domainogram for the supplied scales.
#' @param input A two-column data.frame where the first column represents the fragend and the second the number of UMIs.
#' @param bait_coordintes GRanges object with the coordinates of the bait used for the experiment.
#' @param scales Scales to use for smoothing and plotting the domainogram. Default: 1:150.
#' @param bait_exclusion Region around the bait (in bp) to be excluded from the analysis. Default: 3kb.
#' @param bait_upstream Number of bp upstream of the bait to use for the analysis. Default: 500kb.
#' @param bait_downstream Number of bp downstream of the bait to use for the analysis. Default: 500kb.
#' @export
createDomainogram <- function(input,
                              bait_coordinates,
                              scales=5:150,
                              bait_exclusion=3e3,
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

  ## Divide upstream and downstream coordintes
  upstream <- input[,1] < GenomicRanges::start(bait_coordinates)
  upstream <- seq(1, nrow(input))[upstream]
  downstream <- input[,1] > GenomicRanges::start(bait_coordinates)
  downstream <- seq(1, nrow(input))[downstream]

  ## Get cumsum separately upstream and downstream
  cumsum_up <- cumsum(input[upstream,2])
  cumsum_down <- cumsum(input[downstream,2])

  ## Create domainogram
  dgram <- matrix(input[,1])
  for (s in scales) {
    dgram_up <- c(rep(NA, s),
                  (cumsum_up[(2*s+1):upstream[length(upstream)]] -
                     cumsum_up[1:(upstream[length(upstream)-2*s])])/(2*s),
                  rep(NA,s))

    dgram_dw <- c(rep(NA, s),
                  (cumsum_down[(2*s+1):length(downstream)] -
                     cumsum_down[1:(length(downstream)-2*s)])/(2*s),
                  rep(NA,s))

    dgram <- cbind(dgram, c(dgram_up, dgram_dw))
  }

  colnames(dgram) <- c("start", scales)

  return(dgram)
}
