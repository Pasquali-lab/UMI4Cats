#' Create Domainogram
#'
#' Using as input the raw UMIs, this function creates a domainogram for the
#' supplied scales.
#' @param umi4c \linkS4class{UMI4C} object as generated by \code{\link{makeUMI4C}}.
#' @param scales Integer vector indicating the number of scales to use for the
#' domainogram creation. Default: 5:150.
#' @param normalized Logical whether the the resulting domainograms should be
#' normalized or not. Default: TRUE.
#' @return A matrix where the first column represents the fragment end
#' coordinates (start) and the rest represent the number of UMIs found when
#' using a specific scale.
calculateDomainogram <- function(umi4c,
    scales = 5:150,
    normalized = TRUE) {

    ## Get cumsum separately upstream and downstream
    cumsum_up <- apply(assay(umi4c)[rowRanges(umi4c)$position == "upstream", , drop = FALSE], 2, cumsum)
    cumsum_down <- apply(assay(umi4c)[rowRanges(umi4c)$position == "downstream", , drop = FALSE], 2, cumsum)

    ## Create domainogram
    for (i in seq_len(ncol(assay(umi4c)))) {
        dgram_up <- matrix(NA, nrow = nrow(cumsum_up), ncol = length(scales))
        dgram_dw <- matrix(NA, nrow = nrow(cumsum_down), ncol = length(scales))

        min_scale <- scales[1]
        for (s in scales) {
            dgram_up[, s - (min_scale - 1)] <- c(
                rep(NA, s),
                (cumsum_up[(2 * s + 1):nrow(cumsum_up), i] -
                    cumsum_up[seq_len(nrow(cumsum_up) - 2 * s), i]) / (2 * s),
                rep(NA, s)
            )

            dgram_dw[, s - (min_scale - 1)] <- c(
                rep(NA, s),
                (cumsum_down[(2 * s + 1):nrow(cumsum_down), i] -
                    cumsum_down[seq_len(nrow(cumsum_down) - 2 * s), i]) / (2 * s),
                rep(NA, s)
            )
        }

        dgram <- rbind(
            dgram_up,
            dgram_dw
        )
        rownames(dgram) <- rownames(assay(umi4c))
        colnames(dgram) <- scales

        dgram(umi4c)[[colnames(assay(umi4c))[i]]] <- dgram
    }

    ## Normalize dgrams
    if (normalized) {
        dgram_norm <- dgram(umi4c)

        for (i in seq_len(length(dgram_norm))) {
            dgram_norm[[i]] <- dgram_norm[[i]] * assays(umi4c)$norm_mat[, i]
        }

        dgram(umi4c) <- dgram_norm
    }

    return(umi4c)
}

#' Get normalization matrix
#'
#' Will return a normalization matrix.
#' @param umi4c \linkS4class{UMI4C} object as generated by \code{\link{makeUMI4C}}.
#' @param norm_bins Numeric vector with the genomic bins to use for
#' normalization. Default: 1K, 10K, 100K, 1Mb.
#' @param post_smooth_win Numeric indicating the smoothing window to use.
#' Default: 50.
#' @param r_expand Numeric indicanting the expansion value for normalization.
#' Default: 1.2.
#' @return Creates a matrix of normalization factors using as a reference the
#' profile specified in the \linkS4class{UMI4C} object.
getNormalizationMatrix <- function(umi4c,
    norm_bins = 10^(3:6),
    post_smooth_win = 50,
    r_expand = 1.2) {
    bins <- sort(norm_bins)

    bait_start <- GenomicRanges::start(bait(umi4c))
    coords <- start(rowRanges(umi4c))
    dist <- abs(coords - bait_start) # Ge distances to bait

    # Create empty matrix for normalization
    norm_mat <- matrix(NA,
        ncol = ncol(assay(umi4c)),
        nrow = nrow(assay(umi4c))
    )

    for (s in seq_len(ncol(assay(umi4c)))) { # For each sample present
        dgram <- metadata(umi4c)$dgram[[s]] # get domainogram
        mols <- assay(umi4c)[, s] # get num UMIs

        # Create vector containing total UMIs in each fragment
        # classified by their distance to bait
        sum_vec <- rep(NA, length(coords))

        # First bin, 1kb
        f_bin <- dist < bins[1]
        tot_in_bin <- sum(mols[f_bin]) # stat=linear
        sum_vec[f_bin] <- tot_in_bin

        # Rest of the bins
        for (i in 2:length(bins)) {
            f_bin_exp <- dist < bins[i] * r_expand & dist >= bins[i - 1] / r_expand
            f_bin <- dist < bins[i] & dist >= bins[i - 1]

            tot_in_bin <- sum(mols[f_bin_exp])
            sum_vec[f_bin] <- tot_in_bin
        }

        # Distances > than last bin
        f_bin <- dist >= bins[length(bins)]
        tot_in_bin <- sum(mols[f_bin])
        sum_vec[f_bin] <- tot_in_bin

        # Linear smoothing of values
        sum_vec <- zoo::rollmean(sum_vec, # Vector of UMIs
            post_smooth_win, # Rolling window width (def 50)
            fill = c(sum_vec[1], NA, sum_vec[length(sum_vec)])
        ) # Filling values

        norm_mat[, s] <- sum_vec # Add vector to general matrix
    }

    colnames(norm_mat) <- colnames(assay(umi4c))
    rownames(norm_mat) <- rownames(assay(umi4c))

    ## Normalize to reference profile
    norm_vector <- norm_mat[, metadata(umi4c)$ref_umi4c]
    norm_mat <- norm_vector * (1 / norm_mat)

    assays(umi4c)$norm_mat <- norm_mat

    return(umi4c)
}

#' Adaptative smoothing of normalized trend
#'
#' Will perform adaptative smoothing will scaling one profile to the reference
#' UMI-4C profile.
#' @param umi4c \linkS4class{UMI4C} object as generated by \code{\link{makeUMI4C}}.
#' @inheritParams makeUMI4C
#' @return Calculates the adaptative trend considering the minimum number of
#' molecules to use for merging different restriction fragments. It also
#' calculates the geometric mean of the coordinates of the merged restriction
#' fragments.
calculateAdaptativeTrend <- function(umi4c,
    sd = 2,
    normalized = TRUE) {
    ## Select min_win_cov
    mols <- colSums(assay(umi4c))
    min_mols <- mols * metadata(umi4c)$min_win_factor

    if (normalized) min_mols[min_mols != min(min_mols)] <- min(min_mols)

    ## Get domainograms and convert NA to 0
    dgrams <- dgram(umi4c)
    dgrams <- lapply(
        dgrams,
        function(x) {
            x[is.na(x)] <- 0
            x
        }
    )

    # Initialize values
    vector_trend <- matrix(NA,
        nrow = nrow(dgrams[[1]]),
        ncol = length(dgrams)
    )

    if (normalized) {
        scale <- rep(NA, nrow(dgrams[[1]]))
        coord <- rep(NA, nrow(dgrams[[1]]))
    } else {
        scale <- matrix(NA,
            nrow = nrow(dgrams[[1]]),
            ncol = length(dgrams)
        )
        coord <- matrix(NA,
            nrow = nrow(dgrams[[1]]),
            ncol = length(dgrams)
        )
    }

    base_coord <- start(rowRanges(umi4c))
    mean_coord <- numeric()

    for (i in seq_len(ncol(dgrams[[1]]))) { # Loop for each scale in dgram
        current_scale <- metadata(umi4c)$scales[i]

        # Get geometric mean coordinates
        mean_coords <- geoMeanCoordinates(
            coords = base_coord,
            scale = current_scale * 2,
            bait_start = GenomicRanges::start(bait(umi4c))
        )

        # Obtain rows in which UMIs > threshold for each sample
        pass <- mapply(function(x, min_win_cov) x[, i] * current_scale * 2 > min_win_cov,
            x = dgrams,
            min_win_cov = min_mols
        )

        if (normalized) {
            # Select regions that are still NA in trend and pass cov threshold on
            # all samples
            sums <- rowSums(pass)
            sel <- rowSums(is.na(vector_trend)) == length(dgrams) & (sums == length(dgrams))

            vector_trend[sel, ] <- do.call(
                cbind,
                lapply(
                    dgrams,
                    function(x) x[sel, i]
                )
            )
            coord[sel] <- mean_coords[sel]
            scale[sel] <- current_scale
        } else {
            tmp_dgram <- do.call(cbind, lapply(dgrams, function(x) x[, i]))

            # Select regions that are still NA in trend and pass cov threshold,
            # sample-wise
            na <- is.na(vector_trend)
            pass <- (na + pass) == 2

            vector_trend[pass] <- tmp_dgram[pass]

            mean_coords_mat <- matrix(rep(mean_coords, length(dgrams)),
                ncol = length(dgrams)
            )
            coord[pass] <- mean_coords_mat[pass]
            scale[pass] <- current_scale
        }
    }

    if (normalized) {
        coord <- matrix(rep(coord, length(dgrams)), ncol = length(dgrams))
        scale <- matrix(rep(scale, length(dgrams)), ncol = length(dgrams))
    }

    base_coord_mat <- matrix(rep(base_coord, length(dgrams)), ncol = length(dgrams))
    coord[is.na(coord)] <- base_coord_mat[is.na(coord)]

    ## Add dimnames to coord and scale
    dimnames(vector_trend) <- dimnames(assay(umi4c))
    dimnames(coord) <- dimnames(assay(umi4c))
    dimnames(scale) <- dimnames(assay(umi4c))

    ## Save info in assays
    assays(umi4c)$trend <- vector_trend
    assays(umi4c)$geo_coords <- coord
    assays(umi4c)$scale <- scale

    ## Add SD
    dev <- sd * sqrt(vector_trend / (scale * 2))
    assays(umi4c)$sd <- dev

    return(umi4c)
}

#' Get geometric mean of given coordinates
#'
#' @param coords Vector of integers representing the coordinates from which to
#' obtain the geometric mean.
#' @param scale Vector of scales indicating how many fragment where merged.
#' @param bait_start Integer indicating the coordinates for the bait start.
#' @return Calculates geometric mean of the provided coordinates, taking into
#' account the distance to
#' the viewpoint and how many restriction fragments are being merged.
geoMeanCoordinates <- function(coords,
    scale,
    bait_start) {
    bait_idx <- which(abs(coords - bait_start) == min(abs(coords - bait_start)))[1]
    offsets <- abs(coords - coords[bait_idx])
    offsets[bait_idx] <- 1 # Avoid 0s in log

    mean_offsets_up <- exp(zoo::rollmean(log(offsets[seq_len(bait_idx)]),
        scale,
        fill = NA
    ))
    mean_offsets_dw <- exp(zoo::rollmean(log(offsets[seq(bait_idx + 1, length(offsets))]),
        scale,
        fill = NA
    ))

    ## Deal with NAs near the bait
    near_bait_up <- vapply(
        seq((bait_idx - scale / 2 + 1), bait_idx),
        function(x) exp(mean(log(offsets[seq(x, bait_idx)]))),
        FUN.VALUE = numeric(1)
    )
    near_bait_down <- vapply(
        seq((bait_idx + 1), (bait_idx + scale / 2 - 1)),
        function(x) exp(mean(log(offsets[seq((bait_idx + 1), x)]))),
        FUN.VALUE = numeric(1)
    )

    mean_offsets_up[(bait_idx - scale / 2 + 1):bait_idx] <- near_bait_up
    mean_offsets_dw[seq_len(scale / 2 - 1)] <- near_bait_down

    mean_offsets_up <- mean_offsets_up * -1
    mean_offsets_up[bait_idx] <- 0

    mean_offsets <- c(mean_offsets_up, mean_offsets_dw)
    mean_coords <- coords[bait_idx] + mean_offsets
    return(mean_coords)
}

