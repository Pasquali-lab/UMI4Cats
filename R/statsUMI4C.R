#' Statistics UMI4C
#'
#' Creates a stats summary file and generates a summary plot describing
#' statistics for processed UMI-4C samples.
#' @inheritParams contactsUMI4C
#' @return Returns a plot summarizing the main statistics of the processed
#' UMI-4C experiments found in \code{wk_dir}. Also, writes a file named
#' \code{stats_summary.txt} in \code{wk_dir/logs} that summarizes all the
#' values represented in the plot.
#' @examples
#' statsUMI4C(wk_dir = system.file("extdata", "CIITA",
#'     package = "UMI4Cats"
#' ))
#' stats <- read.delim(system.file("extdata", "CIITA", "logs", "stats_summary.txt",
#'     package = "UMI4Cats"
#' ))
#' head(stats)
#' @export
statsUMI4C <- function(wk_dir) {
    stats <- createStatsTable(wk_dir)

    utils::write.table(stats[, c(
        "sample_id",
        "specific_reads", "nonspecific_reads",
        "filtered_reads", "filtout_reads",
        "al_mapped", "al_unmapped",
        "umi"
    )],
    file = file.path(wk_dir, "logs", "stats_summary.txt"),
    row.names = FALSE,
    quote = FALSE,
    sep = "\t"
    )

    # Create data.frame for plotting
    counts <- reshape2::melt(stats, stringsAsFactors = FALSE)
    counts$stats <- "Specificity"
    counts$stats[grep("filt", counts$variable)] <- "Filtering"
    counts$stats[grep("al_", counts$variable)] <- "Alignment"
    counts$stats[grep("umi", counts$variable)] <- "UMIs"

    counts$variable <- factor(counts$variable,
        levels = c(
            "specific_reads", "nonspecific_reads",
            "filtered_reads", "filtout_reads",
            "al_mapped", "al_unmapped",
            "umi"
        ),
        labels = c(
            "Specific", "Non-specific",
            "Good quality", "Bad quality",
            "Aligned", "Unaligned",
            "UMIs"
        )
    )
    counts$stats <- factor(counts$stats,
        levels = c(
            "Specificity", "Filtering",
            "Alignment", "UMIs"
        )
    )

    # Generate plot for read stats
    read_stats <-
        ggplot2::ggplot(
            counts[counts$stats != "UMIs", ],
            ggplot2::aes(sample_id, value)
        ) +
        ggplot2::geom_bar(ggplot2::aes(fill = variable),
            stat = "identity",
            position = "fill", color = "black"
        ) +
        ggplot2::scale_fill_manual(
            values = c(
                Specific = "olivedrab3",
                "Non-specific" = "tomato3",
                `Good quality` = "darkolivegreen4",
                `Bad quality` = "grey80",
                Aligned = "deepskyblue3",
                Unaligned = "grey"
            ),
            name = "Read type",
            guide = ggplot2::guide_legend(nrow = 2)
        ) +
        ggplot2::coord_flip() +
        ggplot2::scale_y_continuous(
            name = "Percentage of reads",
            labels = function(x) x * 100
        ) +
        ggplot2::facet_wrap(~stats) +
        cowplot::theme_cowplot() +
        ggplot2::theme(
            legend.position = "top",
            axis.title.y = ggplot2::element_blank()
        )

    # Generate plot for umi stats
    umi_stats <-
        ggplot2::ggplot(
            counts[counts$stats == "UMIs", ],
            ggplot2::aes(sample_id, value)
        ) +
        ggplot2::geom_bar(ggplot2::aes(fill = sample_id),
            stat = "identity",
            color = "black"
        ) +
        ggplot2::geom_label(ggplot2::aes(label = scales::comma(value)),
            hjust = 1, nudge_y = -(0.025 * max(counts$value[counts$stats == "UMIs"]))
        ) +
        ggplot2::coord_flip() +
        ggplot2::scale_y_continuous(
            name = "Number of UMIs",
            labels = function(x) scales::comma(x)
        ) +
        ggplot2::facet_wrap(~stats) +
        cowplot::theme_cowplot() +
        ggplot2::theme(
            legend.position = "none",
            axis.title.y = ggplot2::element_blank()
        )

    # Grid of plots
    cowplot::plot_grid(read_stats,
        umi_stats,
        ncol = 1,
        rel_heights = c(0.6, 0.4)
    )
}

#' Create stats table
#'
#' Create a statistical summary of the UMI-4C experiments analyzed with
#' \code{contactsUMI4C}.
#' @inheritParams statsUMI4C
#' @return Returns a data.frame summarizing all the different statistics
#' for each sample analyzed in \code{wk_dir}.
createStatsTable <- function(wk_dir) {
    # List files in logs
    log_files <- list.files(file.path(wk_dir, "logs"),
        pattern = ".txt", full.names = TRUE
    )
    # Ignore stats_summary if it was already created
    log_files <- log_files[!grepl("stats_summary.txt", log_files)]

    # Read log files
    logs <- lapply(log_files, utils::read.delim, stringsAsFactors = FALSE)

    # Get alignment stats
    al <- logs[[grep("alignment", basename(log_files))]]
    al$sample_id <- gsub("_R[[:digit:]]", "", al$sample_id)

    al_stats <- al %>%
        dplyr::group_by(.data$sample_id) %>%
        dplyr::summarise(
            al_mapped = sum(.data$al_mapped),
            al_unmapped = sum(.data$al_unmapped)
        )

    # Get filtering stats
    spec_stats <- logs[[which(!grepl("alignment", basename(log_files)))]]

    # Get UMI stats
    umi_files <- list.files(file.path(wk_dir, "count"),
        pattern = ".tsv",
        full.names = TRUE
    )
    umis <- lapply(umi_files, utils::read.delim, stringsAsFactors = FALSE)
    limits <- c(
        unique(umis[[1]]$pos_bait) - 1.5e3,
        unique(umis[[1]]$pos_bait) + 1.5e3
    )
    umi_sum <- sapply(umis, function(x) sum(x[x$pos_contact < limits[1] | x$pos_contact > limits[2], 5]))
    names <- unlist(lapply(
        strsplit(basename(umi_files), "_counts.tsv"),
        function(x) x[1]
    ))
    umi_df <- data.frame(
        "sample_id" = names,
        "umi" = umi_sum,
        stringsAsFactors = FALSE
    )

    # Merge data
    stats <- dplyr::left_join(spec_stats, al_stats)
    stats$nonspecific_reads <- stats$total_reads - stats$specific_reads
    stats$filtout_reads <- stats$specific_reads - stats$filtered_reads
    stats <- stats[, -which(colnames(stats) == "total_reads")]
    stats <- dplyr::left_join(stats, umi_df)

    return(stats)
}

#' Summarize BAM file
#'
#' Get summary of interesting bam statistics
#' @param bam_file Path for the bam file.
#' @param mapped Logical indicating whether to extract mapped reads.
#' @param secondary Logical indicating whether to extract secondary aligned
#' reads.
#' @return Returns a numeric containing the number of reads in \code{bam_file}
#' that has the specified \code{mapped} and \code{secondary} status.
.getSummaryBam <- function(bam_file, mapped = TRUE, secondary = FALSE) {
    reads <- Rsamtools::countBam(bam_file,
        param = Rsamtools::ScanBamParam(flag = Rsamtools::scanBamFlag(
            isPaired = FALSE,
            isUnmappedQuery = !mapped,
            isSecondaryAlignment = secondary
        ))
    )$records

    return(reads)
}
