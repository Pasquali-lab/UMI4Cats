#' Statistics UMI4C
#'
#' Creates a stats file and generates a summary plot describing interesting statistics
#' for processed UMI-4C samples.
#' @inheritParams contactsUMI4C
#' @examples
#' \dontrun{
#' statsUMI4C(fastq_dir="raw_fastq",
#'            wk_dir="SOCS1")
#' }
#' @export
statsUMI4C <- function(fastq_dir,
                       wk_dir) {
  # Select files necessary for stats
  raw_files <- list.files(fastq_dir,
                          pattern = '_R1.fastq$',
                          full.names=TRUE)
  wk_files <- list.files(file.path(wk_dir, 'prep'),
                         pattern = '_prefiltered_R1.fastq',
                         full.names=TRUE)
  bam <- list.files(file.path(wk_dir, "alignment"),
                    pattern=".bam$",
                    full.names=TRUE)
  samples <- gsub("_R1.fastq", "", basename(raw_files))

  # Generate stats df for filtering and alignment
  counts_df <-
    lapply(samples,
           function(x) data.frame(sample=x,
                                  filt_raw=R.utils::countLines(raw_files[grep(x, raw_files)]),
                                  filt_pass=R.utils::countLines(wk_files[grep(x, wk_files)]),
                                  al_mapped=.getSummaryBam(bam[grep(x, bam)], mapped=TRUE),
                                  al_unmapped=.getSummaryBam(bam[grep(x, bam)], mapped=FALSE),
                                  al_secondary=.getSummaryBam(bam[grep(x, bam)], mapped=TRUE, secondary=TRUE)))
  counts_df <- do.call(rbind, counts_df)
  counts_df$filt_not_pass <- counts_df$filt_raw - counts_df$filt_pass

  # Add UMI stats
  umi_files <- list.files(file.path(wk_dir, "rst"),
                          pattern="10M.tsv",
                          full.names=TRUE)

  counts_df$umi_nums <- sapply(samples,
                               function(x) sum(read.delim(umi_files[grep(x, umi_files)])[,5]))

  # Write output stats table
  write.table(counts_df,
              file=file.path(wk_dir, "rst", "logs.txt"),
              row.names=FALSE,
              quote=FALSE)

  # Remove unnecessary columns
  counts_df <- counts_df[,-grep("filt_raw", colnames(counts_df))]
  counts_df <- counts_df[,-grep("al_secondary", colnames(counts_df))]

  # Create data.frame for plotting
  counts <- reshape2::melt(counts_df)
  counts$stats <- "Filtering"
  counts$stats[grep("al_", counts$variable)] <- "Alignment"
  counts$stats[grep("umi", counts$variable)] <- "UMIs"

  counts$variable <- factor(counts$variable,
                            levels=c("filt_pass", "filt_not_pass",
                                     "al_mapped", "al_unmapped",
                                     "umi_nums"),
                            labels=c("Specific", "Non-specific",
                                     "Aligned", "Unaligned", "UMIs"))
  counts$stats <- factor(counts$stats,
                         levels=c("Filtering", "Alignment", "UMIs"))

  # Generate plot for read stats
  read_stats <-
    ggplot2::ggplot(counts[counts$stats!="UMIs",],
                    ggplot2::aes(sample, value)) +
    ggplot2::geom_bar(ggplot2::aes(fill=variable),
                      stat="identity",
                      position="fill", color="black") +
    ggplot2::scale_fill_manual(values=c(Specific="olivedrab3", "Non-specific"="tomato3",
                                        Aligned="deepskyblue3", Unaligned="grey"),
                               name="Read type",
                               guide=ggplot2::guide_legend(nrow=1)) +
    ggplot2::coord_flip() +
    ggplot2::scale_y_continuous(name="Percentage of reads",
                                labels=function(x) x*100) +
    ggplot2::facet_wrap(~stats) +
    cowplot::theme_cowplot() +
    ggplot2::theme(legend.position="top",
                   axis.title.y=ggplot2::element_blank())

  # Generate plot for umi stats
  umi_stats <-
    ggplot2::ggplot(counts[counts$stats=="UMIs",],
                    ggplot2::aes(sample, value)) +
    ggplot2::geom_bar(ggplot2::aes(fill=sample), stat="identity",
                      color="black") +
    ggplot2::geom_label(ggplot2::aes(label=scales::comma(value)),
                        hjust=1, nudge_y=-(0.025*max(counts$value[counts$stats=="UMIs"]))) +
    ggplot2::coord_flip() +
    ggplot2::scale_y_continuous(name="Number of UMIs",
                                labels=function(x) scales::comma(x)) +
    ggplot2::facet_wrap(~stats) +
    cowplot::theme_cowplot() +
    ggplot2::theme(legend.position="none",
                   axis.title.y=ggplot2::element_blank())

  # Grid of plots
  cowplot::plot_grid(read_stats,
                     umi_stats,
                     ncol=1,
                     rel_heights = c(0.6, 0.4))

}

#' Summarize BAM file
#'
#' Get summary of interesting bam statistics
#' @param bam_file Path for the bam file.
#' @param mapped Logical indicating whether to extract mapped reads.
#' @param secondary Logical indicating whether to extract secondary aligned reads.
.getSummaryBam <- function(bam_file, mapped=TRUE, secondary=FALSE) {
  reads <- Rsamtools::countBam(bam_file,
                               param=Rsamtools::ScanBamParam(flag=Rsamtools::scanBamFlag(isPaired=FALSE,
                                                                                         isUnmappedQuery=!mapped,
                                                                                         isSecondaryAlignment=secondary)))$records

  return(reads)
}
