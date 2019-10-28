geoMeanCoordinates <- function(coords,
                               scale,
                               bait_start) {
  bait_idx <- which(abs(coords - bait_start) == min(abs(coords - bait_start)))[1]
  offsets <- abs(coords - coords[bait_idx])
  offsets[bait_idx] <- 1 # Avoid 0s in log

  mean_offsets_up <- exp(zoo::rollmean(log(offsets[1:bait_idx]), scale, fill=NA))
  mean_offsets_dw <- exp(zoo::rollmean(log(offsets[(bait_idx+1):length(offsets)]), scale, fill=NA))

  ## Deal with NAs near the bait
  near_bait_up <- sapply((bait_idx - scale/2 + 1):bait_idx,
                         function(x) exp(mean(log(offsets[x:bait_idx]))))
  near_bait_down <- sapply((bait_idx + 1):(bait_idx + scale/2 - 1),
                           function(x) exp(mean(log(offsets[(bait_idx + 1):x]))))

  mean_offsets_up[(bait_idx - scale/2 + 1):bait_idx] <- near_bait_up
  mean_offsets_dw[1:(scale/2 - 1)] <- near_bait_down

  mean_offsets_up <- mean_offsets_up * -1
  mean_offsets_up[bait_idx] <- 0

  mean_offsets <- c(mean_offsets_up, mean_offsets_dw)
  mean_coords <- coords[bait_idx] + mean_offsets
  return(mean_coords)
}
