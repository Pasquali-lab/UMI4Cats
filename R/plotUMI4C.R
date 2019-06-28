#' Plot a single UMI-4C profile
#'
#' Using a UMI-4C profile object returned from \code{\link{processUMI4C}}, it draws a plot of the profile and its domainogram.
#' @param umi4c_obj Object outputed by the function \code{\link{processUMI4C}} representing a single UMI-4C processed profile.
#' @param ymax Numeric representing the maximum number of UMIs to show in the y axis of the trend plot. Default: Maximum number of UMIs.
#' @param xlim Numeric vector with the upstream and downstream limits of the region to be plotted. Default: Maximum range in the trend object.
#' @param plot_raw Logical indicating wheter to plot the raw UMI-4C counts as points.
#' @param plot_domainogram Logical for plotting the domainogram (not implemented yet).
#' @param protein_coding Show only protein coding genes in the gene track. Default: TRUE.
#' @param highlight If the input is a gene or a list of gene names will highlight the coordinates of such genes.
#' If the input is a GRanges object, will highlight the coordinates represented there.
#' @param highlight_factor Expansion factor for the coordinates of the highlight. The default is 1, meaning that will plot
#' the coordinates/gene as given. When increased, will expand the area of the highlight.
#' @param highlight_alpha Numeric character from 0 to 1 representing the transparency of the highlight. Default: 0.5
#' @export
#' @return A ggplot2 object with the full plot for the UMI-4C profile.
plotUMI4C <- function(umi4c_obj,
                            ymax=NULL,
                            xlim=NULL,
                            plot_raw=T,
                            plot_domainogram=T,
                            protein_coding=T,
                            highlight=NULL,
                            highlight_factor=1,
                            highlight_alpha=0.5) {
  ## TODO: Add parameter to specify rel_heights in plot_grid (if two many genes first row is crowded)

  ## Update parameters
  if (is.null(ymax)) ymax <- ceiling(max(umi4c_obj$trend$trend, na.rm=T))
  if(is.null(xlim)) {
    xlim <- GenomicRanges::GRanges(paste0(as.character(GenomicRanges::seqnames(umi4c_obj$bait)),
                                                          ":", min(umi4c_obj$trend$coord, na.rm=T),
                                                          "-", max(umi4c_obj$trend$coord, na.rm=T)))
    xlim <- c(GenomicRanges::start(xlim),
              GenomicRanges::end(xlim))
  }

  ## Part 1: Trend and raw UMIs ---------
  trendPlot <- plotTrend(umi4c_obj=umi4c_obj,
                         ymax=ymax,
                         xlim=xlim,
                         plot_raw=plot_raw)

  ## Part 2: Domainogram ------------------
  dgramPlot <- plotDgram(umi4c_obj=umi4c_obj,
                         xlim=xlim)

  ## Part 3: Plot gene annotation -------------------
  win <- umi4c_obj$bait
  GenomicRanges::start(win) <- xlim[1]
  GenomicRanges::end(win) <- xlim[2]

  genesPlot <- plotGenes(window=win,
                         protein_coding=protein_coding)

  ## Highlight ------------------------------------
  if (is.character(highlight)) {
    sel <- unique(genesPlot$data[genesPlot$data$external_gene_name==highlight,1:3])
    highlight <- regioneR::toGRanges(sel)
  }

  if(class(highlight)=="GRanges") {
    message("Adding highlight")
    highlight <- GenomicRanges::resize(highlight,
                                       fix="center",
                                       width=GenomicRanges::width(highlight)*highlight_factor)

    high.plot <- annotate("rect",
                                xmin=GenomicRanges::start(highlight),
                                xmax=GenomicRanges::end(highlight),
                                ymin=-Inf,
                                ymax=Inf,
                                fill="goldenrod3", alpha=highlight_alpha, color=NA)

    ## Add highlight to plot
    genesPlot <- genesPlot + high.plot
    trendPlot <- trendPlot + high.plot
  }


  ## Final: Compose the plot ------------------------

  leg.dom <- cowplot::get_legend(dgramPlot)


  main_plot <-
    cowplot::plot_grid(genesPlot + themeXblank(plot.margin=ggplot2::margin(0.5,0,0.5,0, "cm")),
                       trendPlot + themeXblank(plot.margin=ggplot2::margin(0,0,0,0, "cm")),
                       dgramPlot + themeYblank(legend.position="none",
                                                      plot.margin=ggplot2::margin(0,0,0.5,0, "cm")),
                       ncol=1,
                       align="v",
                       rel_heights=c(0.2,0.4,0.4))#c(0.35,0.35,0.3))

  final_plot <- cowplot::plot_grid(main_plot,
                                   leg.dom,
                                   ncol=2,
                                   rel_widths=c(0.75, 0.25))

  return(final_plot)

}

#' Plot single trend
#'
#' @inheritParams plotUMI4C
#' @export
plotTrend <- function(umi4c_obj=umi4c_obj,
                      ymax=NULL,
                      xlim=NULL,
                      plot_raw=T) {
  ## Update parameters
  if (is.null(ymax)) ymax <- ceiling(max(umi4c_obj$trend$trend, na.rm=T))
  if(is.null(xlim)) {
    xlim <- GenomicRanges::GRanges(paste0(as.character(GenomicRanges::seqnames(umi4c_obj$bait)),
                                          ":", min(umi4c_obj$trend$coord, na.rm=T),
                                          "-", max(umi4c_obj$trend$coord, na.rm=T)))
    xlim <- c(GenomicRanges::start(xlim),
              GenomicRanges::end(xlim))
  }

  ymax.exp <- ymax + ymax/10

  ## Plot trend
  trendPlot <-
    ggplot2::ggplot(umi4c_obj$trend,
                    ggplot2::aes(coord, trend)) +
    ggplot2::geom_point(data=umi4c_obj$raw,
                        ggplot2::aes(coord, scales::rescale(UMIs, to=c(0,
                                                                       ceiling(max(umi4c_obj$trend$trend, na.rm=T))))),
                        color="grey", alpha=0.5) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin=devM, ymax=devP),
                         color=NA, fill="dark grey", alpha=0.8) +
    ggplot2::geom_line(ggplot2::aes(group=group)) +
    ggplot2::annotate("point", x=GenomicRanges::start(umi4c_obj$bait), y=ymax,
             pch=25, color="black", fill="black", size=4) +
    ggplot2::annotate("text", x=GenomicRanges::start(umi4c_obj$bait), y=ymax+ymax/20,
             label=GenomicRanges::mcols(umi4c_obj$bait)[1,1],
             fontface=3) +
    ggplot2::scale_y_continuous(name="# Normalized UMIs",
                                limits=c(0, ymax.exp),
                                expand=c(0,0),
                                breaks=scales::pretty_breaks()) +
    scaleXCoordinates(chr=as.character(GenomicRanges::seqnames(umi4c_obj$bait)),
                      limits=xlim,
                      breaks=scales::pretty_breaks(n=4),
                      expand=c(0,0))

  ## Remove raw if not asked for
  if(!plot_raw) trendPlot$layers <- trendPlot$layers[-1]

  return(trendPlot)
}

#' Plot single domainogram
#'
#' @inheritParams plotUMI4C
#' @export
plotDgram <- function(umi4c_obj=umi4c_obj,
                      xlim=NULL) {
  ## Update xlim
  if(is.null(xlim)) {
    xlim <- GenomicRanges::GRanges(paste0(as.character(GenomicRanges::seqnames(umi4c_obj$bait)),
                                          ":", min(umi4c_obj$trend$coord, na.rm=T),
                                          "-", max(umi4c_obj$trend$coord, na.rm=T)))
    xlim <- c(GenomicRanges::start(xlim),
              GenomicRanges::end(xlim))
  }

  ## Modify dgram to long format
  dgram <- as.data.frame(umi4c_obj$dgram)
  dgram$end <- (dgram$start[c(2:nrow(dgram), nrow(dgram))] -
                  dgram$start) + dgram$start

  dgram <- dgram[,c(1,ncol(dgram), 2:(ncol(dgram)-1))]

  dgram.l <- reshape2::melt(dgram,
                            id.vars=1:2,
                            value.vars=3:ncol(umi4c_obj$dgram),
                            variable.name="scales",
                            value.name="value")
  dgram.l$scales <- as.numeric(dgram.l$scales)

  ## Log2 + normalize by maximum value in view
  norm_fact <- max(dgram.l$value, na.rm=T)
  dgram.l$value_norm <- log2(dgram.l$value+1) - log2(norm_fact+1)

  ## Plot dgram
  domainogram <-
    ggplot2::ggplot(dgram.l) +
    ggplot2::geom_rect(ggplot2::aes(xmin=start, xmax=end,
                                    ymin=rev(scales), ymax=rev(scales)+1,
                                    fill=value_norm)) +
    viridis::scale_fill_viridis(option = "B",
                                na.value=NA,
                                labels=function(x) paste0(round(100*2^x), "%"),
                                name="Contacts/Maximum",
                                breaks=scales::pretty_breaks(n=4),
                                guide = ggplot2::guide_colorbar(direction = "horizontal",
                                                                title.position="top",
                                                                barwidth=8)) +
    scaleXCoordinates(chr=as.character(GenomicRanges::seqnames(umi4c_obj$bait)),
                      limits=xlim,
                      breaks=scales::pretty_breaks(n=4),
                      expand=c(0,0))

  return(domainogram)
}
