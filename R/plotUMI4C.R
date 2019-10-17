plotUMI4C <- function(umi4c,
                      grouping="condition",
                      dgram_function="quotient",
                      dgram_plot=TRUE,
                      colors=NULL,
                      xlim=NULL,
                      ylim=NULL,
                      protein_coding=TRUE,
                      longest=TRUE,
                      rel_heights=c(0.25, 0.4, 0.25),
                      font_size=14) {

  if (is.null(xlim)) {
    xlim <- c(GenomicRanges::start(metadata(umi4c)$region),
            GenomicRanges::end(metadata(umi4c)$region))
  }

    ## Get colors
  factors <- unique(colData(umi4c)[,grouping])
  if (class(factors)=="DataFrame") factors <- do.call(paste, colData(umi4c)[,grouping])

  if (is.null(colors)) colors <- getColors(factors)

  if (length(dgram(umi4c))==1 | length(factors)>2) dgram_plot <- FALSE

  trend_plot <- plotTrend(umi4c,
                          grouping=grouping,
                          xlim=xlim,
                          ylim=ylim,
                          colors=colors)


  genes_plot <- plotGenes(window=metadata(umi4c)$region,
                          protein_coding=protein_coding,
                          longest=longest,
                          xlim=xlim)

  if (dgram_plot) {
    dgram_plot <- plotDomainogram(umi4c,
                                  grouping=grouping,
                                  dgram_function=dgram_function,
                                  colors=colors,
                                  xlim=xlim) + cowplot::theme_cowplot(font_size)

    trend_theme <- cowplot::theme_cowplot(font_size) + themeXblank(legend.position="bottom",
                                                          legend.justification="center")
  } else {
    dgram_plot <- NULL
    trend_theme <- cowplot::theme_cowplot(font_size) + ggplot2::theme(legend.position="bottom",
                                                             legend.justification="center")
  }

  umi4c_plot <- list(genes_plot + cowplot::theme_nothing(font_size),
                     trend_plot + trend_theme,
                     dgram_plot + ggplot2::theme(legend.position="bottom"))

  ## Remove dgram data if dgram_plot is false
  umi4c_keep <- !sapply(umi4c_plot, is.null)

  if (any(!umi4c_keep)) {
    umi4c_plot <- umi4c_plot[umi4c_keep]
    rel_heights <- rel_heights[umi4c_keep]
  }

  ## Extract legends and plot them separately
  legends <- lapply(umi4c_plot[-1], cowplot::get_legend)
  legends_plot <- cowplot::plot_grid(plotlist=legends, nrow=1, align="h")

  ## Remove legends from plot
  umi4c_plot <- lapply(umi4c_plot, function(x) x + ggplot2::theme(legend.position="none"))

  ## Plot main
  main_plot <- cowplot::plot_grid(plotlist=umi4c_plot,
                                  ncol=1,
                                  align="v",
                                  rel_heights=rel_heights)



  cowplot::plot_grid(legends_plot, main_plot,
                     ncol=1,
                     rel_heights = c(0.15,0.85))

}

plotDomainogram <- function(umi4c,
                            grouping="condition",
                            dgram_function="quotient", # or "difference"
                            colors,
                            xlim=NULL) {
  factor <- unique(colData(umi4c)[, grouping])
  if (class(factor)=="DataFrame") factor <- do.call(paste, colData(umi4c)[,grouping])


  if (length(factor)>2) stop("Error in 'plotDomainogram':\n
                             dgram_grouping' cannot have more than two levels. Choose another
                             variable for grouping or refactor the column to only have two levels.")

  ## TODO: Consider length(factor)==1: plot single domainogram
  dgram <- dgram(umi4c)

  ## Sum dgrams from same factor
  ids_1 <- colData(umi4c)$sampleID[grep(factor[1], colData(umi4c)[,grouping])]
  ids_2 <- colData(umi4c)$sampleID[grep(factor[2], colData(umi4c)[,grouping])]

  dgram_merged <- list()
  dgram_merged[[factor[1]]] <- Reduce('+', dgram[ids_1])
  dgram_merged[[factor[1]]][is.na(dgram_merged[[factor[1]]])] <- 0

  dgram_merged[[factor[2]]] <- Reduce('+', dgram[ids_2])
  dgram_merged[[factor[2]]][is.na(dgram_merged[[factor[2]]])] <- 0

  ## Create dgram of difference
  if (dgram_function=="difference") {
    dgram_diff <- log2(1 + dgram_merged[[factor[2]]]) - log2(1 + dgram_merged[[factor[1]]])
    lab_legend <- "Diff "
  } else if (dgram_function=="quotient") {
    dgram_diff <- log2(dgram_merged[[factor[2]]]/dgram_merged[[factor[1]]])
    lab_legend <- "Quot "
  }

  ## Create melted dgram
  dgram_diff <- reshape2::melt(dgram_diff)
  colnames(dgram_diff) <- c("contact_id", "scales", "value")

  ## Add coordinates
  dgram_diff$start <- rep(GenomicRanges::start(umi4c),
                          length(unique(dgram_diff$scales)))
  dgram_diff$end <- rep((GenomicRanges::start(umi4c)[c(2:length(umi4c), length(umi4c))] -
                           GenomicRanges::start(umi4c)) + GenomicRanges::start(umi4c),
                        length(unique(dgram_diff$scales)))


  dgram_plot <-
    ggplot2::ggplot(dgram_diff) +
    ggplot2::geom_rect(ggplot2::aes(xmin=start, xmax=end,
                                    ymin=scales, ymax=scales+1,
                                    fill=value)) +
    ggplot2::scale_fill_gradientn(colors=c(darken(colors[1], factor=10),
                                           colors[1], "white",
                                           colors[2],
                                           darken(colors[2], factor=10)),
                                  na.value = NA,
                                  name=as.expression(bquote("Dgram "*.(lab_legend)*log[2]*" UMIs")),
                                  breaks=scales::pretty_breaks(n=4),
                                  guide = ggplot2::guide_colorbar(direction = "horizontal",
                                                                  title.position="top",
                                                                  barwidth=8)) +
    ggplot2::scale_y_reverse(name="",
                             breaks=c(min(metadata(umi4c)$scales),
                                     max(metadata(umi4c)$scales)),
                             expand=c(0,0)) +
    ggplot2::scale_x_continuous(labels=function(x) round(x/1e6,2),
                                name=paste0("Coordinates",
                                            GenomicRanges::seqnames(metadata(umi4c)$bait),
                                            "(Mb)")) +
    ggplot2::coord_cartesian(xlim=xlim) +
    ggplot2::guides(fill=ggplot2::guide_colorbar(title.position="left",
                                                 label.position="bottom",
                                                 title.vjust=1,
                                                 direction="horizontal"))

  return(dgram_plot)
}


plotTrend <- function(umi4c,
                      grouping="condition",
                      colors,
                      xlim=NULL,
                      ylim=NULL) {
  ## Construct trend df using geo_coords and trend
  trend_df <- data.frame(geo_coord = as.vector(assays(umi4c)$geo_coord),
                         trend = as.vector(assays(umi4c)$trend),
                         sd = as.vector(assays(umi4c)$sd),
                         scale = as.vector(assays(umi4c)$scale),
                          id_contact = rep(rowRanges(umi4c)$id_contact,
                                          ncol(umi4c)))
  trend_df$sample <- rep(colnames(assay(umi4c)),
                         each=nrow(umi4c))

  trend_df <- trend_df[!is.na(trend_df$trend),]

  trend_df <- dplyr::left_join(trend_df,
                               data.frame(colData(umi4c)),
                               by=c(sample="sampleID"))

  if (length(grouping)==1) {
    trend_df <-
      trend_df %>%
      group_by_at(c(grouping, "id_contact")) %>%
      summarise(geo_coord=mean(geo_coord),
                trend=sum(trend),
                sd=mean(sd),
                scale=mean(scale))
  }

  trend_df$relative_position <- "upstream"
  trend_df$relative_position[trend_df$geo_coord>GenomicRanges::start(metadata(umi4c)$bait)] <- "downstream"
  trend_df$grouping_var <- do.call(paste, trend_df[,grouping])

  cols <- rev(legend)
  names(cols) <- NULL

  trend_plot <-
    ggplot2::ggplot(trend_df) +
    ggplot2::geom_ribbon(ggplot2::aes(geo_coord, ymin=trend-sd, ymax=trend+sd,
                                      group=interaction(grouping_var, relative_position),
                                      fill=grouping_var),
                         alpha=0.3, color=NA) +
    ggplot2::geom_line(ggplot2::aes(geo_coord, trend,
                                    group=interaction(grouping_var, relative_position),
                                    color=grouping_var)) +
    ggplot2::scale_color_manual(values=colors,
                                name="Trend group") +
    ggplot2::scale_fill_manual(values=colors,
                                name="Trend group") +
    ggplot2::annotate("point", x=GenomicRanges::start(metadata(umi4c)$bait), y=max(ylim),
                      color="black", fill="black", pch=25, size=4) +
    ggplot2::coord_cartesian(xlim=xlim, ylim=ylim) +
    ggplot2::scale_y_continuous(name="UMIs normalized trend",
                                breaks=scales::pretty_breaks(),
                                expand=c(0,0)) +
    ggplot2::scale_x_continuous(labels=function(x) round(x/1e6,2),
                                name=paste0("Coordinates",
                                            GenomicRanges::seqnames(metadata(umi4c)$bait),
                                            "(Mb)")) +
    ggplot2::theme(legend.position="bottom")

  return(trend_plot)
}

#' Plot genes
#'
#' @param window GRanges object with coordinates from which to plot existing genes.
#' @param protein_coding Plot only protein coding genes. Default: TRUE
#' @export
#' @importFrom magrittr %>%
plotGenes <- function(window,
                      protein_coding=TRUE,
                      longest=TRUE,
                      xlim=NULL) {
  ## Get gene names in region
  genes_sel <- unique(IRanges::subsetByOverlaps(hg19_gene_annoation_ensemblv75, window))

  if(protein_coding) genes_sel <- genes_sel[genes_sel$gene_biotype=="protein_coding",]
  if(longest) {
    tx_ids_sel <- unique(genes_sel$tx_id[genes_sel$longest])
    genes_sel <- genes_sel[genes_sel$tx_id %in% tx_ids_sel]
  }
  ## Edit genes
  distance <- GenomicRanges::width(window)*0.01

  ## Add stepping
  genes_step <- addStepping(genes_sel[genes_sel$type=="GENE",], window, 2)
  genes_uni <- data.frame(genes_step)

  genes_exon <- data.frame(genes_sel[genes_sel$type=="EXON",])
  genes_exon <- dplyr::left_join(genes_exon,
                                 genes_uni[,c(8,12)],
                                 by=c(tx_id="tx_id"))


  ## Plot genes--------------
  genesPlot <-
    ggplot2::ggplot(data=genes_uni) +
    ggplot2::geom_segment(data=genes_uni,
                          ggplot2::aes(x=start, y=stepping,
                                       xend=end, yend=stepping)) +
    ggplot2::geom_rect(data=genes_exon,
                       ggplot2::aes(xmin=start, xmax=end,
                                    ymin=(stepping-0.3), ymax=(stepping+0.3)),
                       fill="grey39", color="grey39") +
    ggplot2::geom_text(data=genes_uni,
                       ggplot2::aes(x=end, y=stepping, label=gene_name),
                       colour="black",
                       hjust=0, fontface=3, nudge_x=distance,
                       size=3) +
    ggplot2::coord_cartesian(xlim=xlim)

  return(genesPlot)

}

#' Add stepping for plotting genes
#'
#' Given a GRanges dataset representing genes, will add an arbitrary value for them to be plotted in
#' the Y axis without overlapping each other.
#' @param genesDat GRanges object containing gene information.
#' @param coordinates GRanges object with coordinates you want to plot.
#' @import GenomicRanges
addStepping <- function(genesDat,
                        coordinates,
                        mcol.name) {
  ## Create extension for avoiding overlap with gene names
  ext <- sapply(mcols(genesDat)[,mcol.name], nchar) * width(coordinates)/30
  genesDat.ext <- regioneR::extendRegions(genesDat, extend.end=ext)

  ## Add stepping to data
  genesDat$stepping <- disjointBins(genesDat.ext,
                                    ignore.strand=TRUE)

  return(genesDat)
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


#' Theme Y blank
#' @export
themeXblank <- function(...) {
  ggplot2::theme(axis.text.x=ggplot2::element_blank(),
                 axis.title.x=ggplot2::element_blank(),
                 axis.line.x=ggplot2::element_blank(),
                 axis.ticks.x=ggplot2::element_blank(),
                 ...)
}


getColors <- function(factors) {
    if (length(factors)==2) {
      colors <- c("darkorchid3", "darkorange3")
      # names(colors) <- factors
    } else if (length(factors)>2) {
      colors <- RColorBrewer::brewer.pal(n=length(factors), name="Set1")
      # names(colors) <- colors
    } else if (length(factors)==1) {
      colors <- "darkorchid3"
      # names(colors) <- factors
    }
  return(colors)
}
