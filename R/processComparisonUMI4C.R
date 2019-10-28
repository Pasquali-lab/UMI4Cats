#' Compare two UMI-4C profiles
#'
#' Scales the umi4c_obj to the ref_umi4c_obj, creating all necessary elements to produce the comparative plot.
#' @param umi4c_obj Object outputed by the function \code{\link{processUMI4C}} representing a single UMI-4C processed profile.
#' @param ref_umi4c_obj Object outputed by the function \code{\link{processUMI4C}} representing a single UMI-4C processed profile to use as reference.
#' \code{umi4c_obj} will be scaled to the \code{ref_umi4c_obj}.
#' @param min_win_cov Minimum number of UMIs requiered in a region to perform the adaptaive smoothing. Default: NULL,
#' number is calculated by the function \code{\link{getMinWinMols}} using the factor specified below.
#' @param norm_bins Bins used to normalize the number of UMIs from \code{umi4c_obj} to \code{umi4c_ref}. Default: 10^(3:6).
#' @param factor Factor representing the proportion of UMIs from the total UMIs in the analyzed region needed as a minimum
#' to perform the adaptative smoothing.
#' @param sd Standard devition value for the adaptative smoothing. Default: 2 (extracted from Schwartzman et al).
#' @export
#' @return A list object with the following named objects:
processComparisonUMI4C <- function(umi4c_obj,
                                   ref_umi4c_obj,
                                   min_win_cov=NULL,
                                   norm_bins=10^(3:6),
                                   factor=0.025,
                                   sd=2) {

  ## Normalize domainogram from treat to ref
  umi4c_obj <- normalizeDomainogram(umi4c_obj,
                                    ref_umi4c_obj,
                                    norm_bins=norm_bins)

  ## Summarize raw
  raw <- summarizeRaw(umi4c_obj,
                      ref_umi4c_obj)

  ## Create domainogram of comparison
  dif_dgram <- getDomainogramComparison(umi4c_obj,
                                        ref_umi4c_obj)

  ## Get default min_win_cov
  if (is.null(min_win_cov)) min_win_cov <- getMinWinMols(sum(ref_umi4c_obj$raw[,2]),
                                                         factor=factor)

  ## Plot comparative trend
  trend_comp <- smoothTrendAdaptativeComp(umi4c_obj,
                                          ref_umi4c_obj,
                                          min_win_cov=min_win_cov,
                                          sd=sd)

  trend_comp$dgram <- dif_dgram
  trend_comp$raw <- raw

  return(trend_comp)
}

#' Summarize raw UMIs
#'
#' @inheritParams processComparisonUMI4C
#' @export
summarizeRaw <- function(umi4c_obj,
                         ref_umi4c_obj) {
  raw <- umi4c_obj$raw
  raw$type <- "Treat"

  raw_ref <- ref_umi4c_obj$raw
  raw_ref$type <- "Ref"

  raw <- rbind(raw,
               raw_ref)

  return(raw)
}

#' Get comparative domainogram
#'
#' Will return values of a comparative domainogram representing the differences between log2 UMIs.
#' @inheritParams processComparisonUMI4C
#' @export
getDomainogramComparison <- function(umi4c_obj,
                                     ref_umi4c_obj) {
  dgramNorm <- umi4c_obj$normDgram$dgram
  # dgramNorm[is.na(dgramNorm)] <- 0
  dgramRef <- ref_umi4c_obj$dgram
  # dgramRef[is.na(dgramRef)] <- 0

  ## Create dgram of difference
  dif_dgram <- log2(1 + dgramNorm[,-1]) - log2(1 + dgramRef[,-1])
  # dif_dgram <- log2(dgramNorm[,-1]/dgramRef[,-1])

  dif_dgram <- cbind(dgramNorm[,1], dif_dgram)
  colnames(dif_dgram)[1] <- "start"

  return(dif_dgram)
}

#' Adaptative smoothing of scaled trend
#'
#' Will perform adaptative smoothing will scaling one profile to the reference UMI-4C profile.
#' @inheritParams processComparisonUMI4C
#' @export
smoothTrendAdaptativeComp <- function(umi4c_obj,
                                      ref_umi4c_obj,
                                      min_win_cov=50,
                                      sd=2) {
  bait_coordinates <- umi4c_obj$bait
  norm_factors <- umi4c_obj$normDgram$factors

  # Get domainograms
  dgram <- umi4c_obj$dgram
  dgram[is.na(dgram)] <- 0
  dgramRef <- ref_umi4c_obj$dgram
  dgramRef[is.na(dgramRef)] <- 0

  # Initialize values
  vector <- rep(NA, nrow(dgram))
  vector_ref <- rep(NA, nrow(dgramRef))
  scale <- rep(NA, nrow(dgram))
  coord <- rep(NA, nrow(dgram))

  base_coord <- dgram[,1]
  mean_coord <- numeric()
  base_scales <- as.numeric(colnames(dgram)[-1])

  for (i in 2:ncol(dgram)) {
    cur_scale <- as.numeric(colnames(dgram)[i])

    f <- is.na(vector) &
      (dgram[,i] * cur_scale * 2 > min_win_cov &
         dgramRef[,i] * cur_scale * 2 > min_win_cov)

    vector[f] <- dgram[f,i]
    vector_ref[f] <- dgramRef[f,i]

    ## Get geometric mean coordinates
    mean_coords <- geoMeanCoordinates(coords=base_coord,
                                      scale=cur_scale*2,
                                      bait_start=GenomicRanges::start(bait_coordinates))
    coord[f] <- mean_coords[f]
    scale[f] <- cur_scale
  }

  vector_norm <- vector * norm_factors
  coord[is.na(coord)] <- base_coord[is.na(coord)]

  ## Create trend TREAT
  trend <- data.frame(coord=coord,
                      trend=vector_norm,
                      scale=scale)
  trend$group <- "upstream"
  trend$group[trend$coord > GenomicRanges::start(bait_coordinates)] <- "downstream"
  trend <- unique(trend)

  ## Add SD
  dev <- sd * sqrt(trend$trend/(trend$scale*2))
  trend$devP <- trend$trend + dev
  trend$devM <- trend$trend - dev
  trend$type <- "Treat"

  ## Create trend REF
  trend_ref <- data.frame(coord=coord,
                          trend=vector_ref,
                          scale=scale)
  trend_ref$group <- "upstream"
  trend_ref$group[trend_ref$coord > GenomicRanges::start(bait_coordinates)] <- "downstream"
  trend_ref <- unique(trend_ref)

  ## Add SD
  dev <- sd * sqrt(trend_ref$trend/(trend_ref$scale*2))
  trend_ref$devP <- trend_ref$trend + dev
  trend_ref$devM <- trend_ref$trend - dev
  trend_ref$type <- "Ref"

  ## Merge trends
  trend <- rbind(trend, trend_ref)



  ## Return comparison object
  comp <- list(trend=trend,
               bait=bait_coordinates)

  return(comp)
}

#' Normalize domainogram to reference
#'
#' Will normalize the domainogram of a given umi4c_obj to the reference.
#' @inheritParams processComparisonUMI4C
#' @export
normalizeDomainogram <- function(umi4c_obj,
                                 ref_umi4c_obj,
                                 norm_bins=10^(3:6)) {

  norm1 <- getNormalizationVector(umi4c_obj,
                                  norm_bins=norm_bins)
  norm2 <- getNormalizationVector(ref_umi4c_obj,
                                  norm_bins=norm_bins)

  norm_vector <- norm2/norm1

  normDgram <- umi4c_obj$dgram
  normDgram[,-1] <- normDgram[,-1] * norm_vector

  umi4c_obj$normDgram <- list(dgram=normDgram,
                              factors=norm_vector)

  return(umi4c_obj)
}

#' Get normalization vector
#'
#' Will return a vector for normalizing domainogram.
#' @inheritParams processComparisonUMI4C
#' @export
getNormalizationVector <- function(umi4c_obj,
                                   norm_bins=10^(3:6),
                                   post_smooth_win=50,
                                   r_expand=1.2) {
  bins <- sort(norm_bins)

  bait_start <- GenomicRanges::start(umi4c_obj$bait)

  dgram <- umi4c_obj$dgram
  mols <- umi4c_obj$raw$UMIs
  coords <- umi4c_obj$dgram[,1]
  dist <- abs(coords - bait_start)

  sum_vec <- rep(NA, length(coords))

  f_bin <- dist < bins[1]
  tot_in_bin <- sum(mols[f_bin]) # stat=linear
  sum_vec[f_bin] <- tot_in_bin

  for (i in 2:length(bins)) {
    f_bin_exp <- dist < bins[i]*r_expand & dist >= bins[i-1]/r_expand
    f_bin <- dist < bins[i] & dist >= bins[i-1]

    tot_in_bin <- sum(mols[f_bin_exp])
    sum_vec[f_bin] <- tot_in_bin

    message("bin ", bins[i], " total ", tot_in_bin)
  }

  f_bin <- dist >= bins[length(bins)]
  tot_in_bin <- sum(mols[f_bin])
  sum_vec[f_bin] <- tot_in_bin

  sum_vec <- zoo::rollmean(sum_vec,
                           post_smooth_win,
                           fill=c(sum_vec[1], NA, sum_vec[length(sum_vec)]))

  return(sum_vec)
}
