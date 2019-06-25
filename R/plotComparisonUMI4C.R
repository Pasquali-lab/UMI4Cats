#' Comparative plot of UMI4C profiles
#'
#' Using a comparison object returned from \code{\link{processTrendComparison}}, it draws a comparative plot of the
#' two profiles and a domainogram of the difference.
#' @param comp_obj Object outputed by the function \code{\link{processTrendComparison}}.
#' @param ymax Numeric representing the maximum number of UMIs to show in the y axis of the trend plot. Default: Maximum number of UMIs.
#' @param xlim Numeric vector with the upstream and downstream limits of the region to be plotted. Default: Maximum range in the trend object.
#' @param plot_domainogram Logical for plotting the domainogram (not implemented yet).
#' @param protein_coding Show only protein coding genes in the gene track. Default: TRUE.
#' @param legend Named character vector with the names to use for the legend and colors to use for plotting. The first
#' element should be the "Treatment" and the second one the sample used as reference.
#' @param highlight If the input is a gene or a list of gene names will highlight the coordinates of such genes.
#' If the input is a GRanges object, will highlight the coordinates represented there.
#' @param highlight_factor Expansion factor for the coordinates of the highlight. The default is 1, meaning that will plot
#' the coordinates/gene as given. When increased, will expand the area of the highlight.
#' @param highlight_alpha Numeric character from 0 to 1 representing the transparency of the highlight. Default: 0.5
#' @export
#' @return A ggplot2 object with the full comparative plot.
plotComparisonUMI4C <- function(comp_obj,
                                ymax=NULL,
                                xlim=NULL,
                                plot_domainogram=T,
                                protein_coding=T,
                                legend=c("Treatment"="darkorange3",
                                         "Reference"="darkorchid3"),
                                highlight=NULL,
                                highlight_factor=1,
                                highlight_alpha=0.5) {
  ## Update parameters
  if (is.null(ymax)) ymax <- ceiling(max(comp_obj$trend$trend, na.rm=T))
  if(is.null(xlim)) {
    xlim <- GenomicRanges::GRanges(paste0(as.character(GenomicRanges::seqnames(comp_obj$bait)),
                                          ":", min(comp_obj$trend$coord, na.rm=T),
                                          "-", max(comp_obj$trend$coord, na.rm=T)))
    xlim <- c(GenomicRanges::start(xlim),
              GenomicRanges::end(xlim))
  }

  ## Part 1: Trend and raw UMIs ---------
  trendPlot <- plotTrendComp(comp_obj=comp_obj,
                             ymax=ymax,
                             xlim=xlim,
                             legend=legend)

  ## Part 2: Domainogram ------------------
  dgramPlot <- plotDgramComp(comp_obj=comp_obj,
                             xlim=xlim,
                             legend=legend)

  ## Part 3: Plot gene annotation -------------------
  win <- comp_obj$bait
  GenomicRanges::start(win) <- xlim[1]
  GenomicRanges::end(win) <- xlim[2]

  genesPlot <- plotGenes(window=win,
                         protein_coding=protein_coding)

  ## Add highlight
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
  leg.trend <- cowplot::get_legend(trendPlot)

  plot_legend <- plot_grid(leg.trend,
                           leg.dom,
                           ncol=1)


  main_plot <-
    cowplot::plot_grid(genesPlot + themeXblank(plot.margin=ggplot2::margin(0.5,0,0.5,0, "cm")),
                       trendPlot + themeXblank(legend.position="none",
                                               plot.margin=ggplot2::margin(0,0,0,0, "cm")),
                       dgramPlot + themeYblank(legend.position="none",
                                               plot.margin=ggplot2::margin(0,0,0.5,0, "cm")),
                       ncol=1,
                       align="v",
                       rel_heights=c(0.2,0.4,0.4))

  final_plot <- cowplot::plot_grid(main_plot,
                                   plot_legend,
                                   ncol=2,
                                   rel_widths=c(0.75, 0.25))

  return(final_plot)
}

#' Plot comparative trend
#' @inheritParams plotComparisonUMI4C
#' @export
plotTrendComp <- function(comp_obj,
                      ymax=NULL,
                      xlim=NULL,
                      legend) {
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

  comp_obj$trend$type[grep("Treat", comp_obj$trend$type)] <- names(legend)[1]
  comp_obj$trend$type[grep("Ref", comp_obj$trend$type)] <- names(legend)[2]

  ## Plot trend
  trendPlot <-
    ggplot2::ggplot(comp_obj$trend,
                    ggplot2::aes(coord, trend)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin=devM, ymax=devP, group=type),
                         color=NA, fill="dark grey", alpha=0.8) +
    ggplot2::geom_line(ggplot2::aes(group=interaction(group, type), color=type)) +
    ggplot2::annotate("point", x=GenomicRanges::start(comp_obj$bait), y=ymax,
             pch=25, color="black", fill="black", size=4) +
    ggplot2::annotate("text", x=GenomicRanges::start(comp_obj$bait), y=ymax+ymax/20,
                      label=GenomicRanges::mcols(comp_obj$bait)[1,1],
                      fontface=3) +
    ggplot2::scale_color_manual(values=legend,
                                name="") +
    ggplot2::scale_y_continuous(name="# Normalized UMIs",
                                limits=c(0, ymax.exp),
                                expand=c(0,0),
                                labels=function(x) scales::comma(x),
                                breaks=scales::pretty_breaks()) +
    scaleXCoordinates(chr=as.character(GenomicRanges::seqnames(comp_obj$bait)),
                      limits=xlim,
                      breaks=scales::pretty_breaks(n=4),
                      expand=c(0,0))

  return(trendPlot)
}

#' Plot comparative domainogram
#'
#' @inheritParams plotComparisonUMI4C
#' @export
plotDgramComp <- function(comp_obj,
                          xlim=NULL,
                          legend) {
  ## Update xlim
  if(is.null(xlim)) {
    xlim <- GenomicRanges::GRanges(paste0(as.character(GenomicRanges::seqnames(umi4c_obj$bait)),
                                          ":", min(umi4c_obj$trend$coord, na.rm=T),
                                          "-", max(umi4c_obj$trend$coord, na.rm=T)))
    xlim <- c(GenomicRanges::start(xlim),
              GenomicRanges::end(xlim))
  }

  ## Modify dgram to long format
  dgram <- as.data.frame(comp_obj$dgram)
  dgram$end <- (dgram$start[c(2:nrow(dgram), nrow(dgram))] -
                  dgram$start) + dgram$start

  dgram <- dgram[,c(1,ncol(dgram), 2:(ncol(dgram)-1))]

  dgram.l <- reshape2::melt(dgram,
                            id.vars=1:2,
                            value.vars=3:ncol(umi4c_obj$dgram),
                            variable.name="scales",
                            value.name="value")
  dgram.l$scales <- as.numeric(dgram.l$scales)

  ## Plot dgram
  domainogram <-
    ggplot2::ggplot(dgram.l) +
    ggplot2::geom_rect(ggplot2::aes(xmin=start, xmax=end,
                                    ymin=rev(scales), ymax=rev(scales)+1,
                                    fill=value)) +
    ggplot2::scale_fill_gradientn(colors=c(darken(legend[2], factor=10),
                                           legend[2], "white",
                                           legend[1],
                                           darken(legend[1], factor=10)),
                                  na.value = NA,
                                  name=expression("Diff "*log[2]*" UMIs"),
                                  breaks=scales::pretty_breaks(n=4),
                                  guide = ggplot2::guide_colorbar(direction = "horizontal",
                                                                  title.position="top",
                                                                  barwidth=8)) +
    scaleXCoordinates(chr=as.character(GenomicRanges::seqnames(comp_obj$bait)),
                      limits=xlim,
                      breaks=scales::pretty_breaks(n=4),
                      expand=c(0,0))

  return(domainogram)
}


#' Darken colors
#'
#' @param color Character containing the name or hex value of a color.
#' @param factor Numeric representing a factor by which darken the specified color.
darken <- function(color, factor=1.4){
  col <- col2rgb(color)
  col <- col/factor
  col <- rgb(t(col), maxColorValue=255)
  col
}
