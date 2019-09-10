#' Get normalization matrix
#'
#' Will return a vector for normalizing domainogram.
#' @param umi4c
#' @param norm_bins
#' @param post_smooth_win
#' @param r_expand
#' @export
getNormalizationMatrix <- function(umi4c,
                                   norm_bins=10^(3:6),
                                   post_smooth_win=50,
                                   r_expand=1.2) {
  bins <- sort(norm_bins)

  bait_start <- GenomicRanges::start(metadata(umi4c)$bait)
  coords <- start(rowRanges(umi4c))
  dist <- abs(coords - bait_start)

  norm_mat <- matrix(NA,
                     ncol=ncol(assay(umi4c)),
                     nrow=nrow(assay(umi4c)))

  for (s in 1:ncol(assay(umi4c))) {
    dgram <- metadata(umi4c)$dgram[[s]]
    mols <- assay(umi4c)[,s]

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

    norm_mat[,s] <- sum_vec
  }

  colnames(norm_mat) <- colnames(assay(umi4c))
  rownames(norm_mat) <- rownames(assay(umi4c))

  ## Normalize to reference profile
  norm_vector <- norm_mat[,metadata(umi4c)$parameters$ref_umi4c]
  norm_mat <- norm_vector*(1/norm_mat)

  assays(umi4c)$norm_mat <- norm_mat

  return(umi4c)
}
