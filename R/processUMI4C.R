#' Process a single UMI4C profile
#'
#' This function processes a single UMI4C profile and generates all data necessary for plotting.
#' @param input A two-column data.frame where the first column represents the fragend and the second the number of UMIs.
#' @param bait_coordintes GRanges object with the coordinates of the bait used for the experiment.
#' @param scales Scales to use for smoothing and plotting the domainogram. Default: 1:150.
#' @param bait_exclusion Region around the bait (in bp) to be excluded from the analysis. Default: 3kb.
#' @param bait_upstream Number of bp upstream of the bait to use for the analysis. Default: 500kb.
#' @param bait_downstream Number of bp downstream of the bait to use for the analysis. Default: 500kb.
#' @param min_win_cov Minimum number of UMIs requiered in a region to perform the adaptaive smoothing. Default: NULL,
#' number is calculated by the function getMinWinMols() using the factor specified below.
#' @param factor Factor representing the proportion of UMIs from the total UMIs in the analyzed region needed as a minimum
#' to perform the adaptative smoothing.
#' @param sd Standard devition value for the adaptative smoothing. Default: 2 (extracted from Schwartzman et al).
#' @export
processUMI4C <- function(input,
                         bait_coordinates,
                         scales=5:150,
                         bait_exclusion=3e3,
                         bait_upstream=5e5,
                         bait_downstream=5e5,
                         min_win_cov=NULL,
                         factor=0.025,
                         sd=2) {

  ## Create domainogram
  dgram <- createDomainogram(input=input,
                             bait_coordinates=bait_coordinates,
                             scales=scales,
                             bait_exclusion=bait_exclusion,
                             bait_upstream=bait_upstream,
                             bait_downstream=bait_downstream)

  ## Get default min_win_cov
  if (is.null(min_win_cov)) min_win_cov <- getMinWinMols(sum(input[,2]),
                                                         factor=factor)

  ## Create adaptative smoothed trend
  trend <- smoothTrendAdaptative(dgram,
                                 min_win_cov=min_win_cov,
                                 bait_coordinates=bait_coordinates,
                                 sd=sd)

  ## Obtain raw
  raw <- filterRaw(input=input,
                   bait_coordinates=bait_coordinates,
                   bait_exclusion=bait_exclusion,
                   bait_upstream=bait_upstream,
                   bait_downstream=bait_downstream)

  ## Create UMI4CObject
  umi4c_obj <- list(raw=raw,
                    bait=bait_coordinates,
                    dgram=dgram,
                    trend=trend)

  return(umi4c_obj)
}


# processComparisonUMI4C <- function()
