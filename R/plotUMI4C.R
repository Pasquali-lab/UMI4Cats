plotUMI4C <- function(umi4c,
                      trend_grouping=c("condition", "replicate"),
                      dgram_grouping="condition",
                      dgram_function="quotient",
                      legend=c("Treatment"="darkorange3",
                               "Reference"="darkorchid3"),
                      xlim=NULL,
                      ylim=NULL,
                      protein_coding=TRUE,
                      longest=TRUE) {

  if (is.null(xlim)) {
    xlim <- c(GenomicRanges::start(metadata(umi4c)$region),
            GenomicRanges::end(metadata(umi4c)$region))
  }

  trend_plot <- plotTrend(umi4c,
                          grouping=trend_grouping,
                          xlim=xlim,
                          ylim=ylim)


  genes_plot <- plotGenes(window=metadata(umi4c)$region,
                          protein_coding=protein_coding,
                          longest=longest,
                          xlim=xlim)

  dgram_plot <- plotDomainogram(umi4c,
                                grouping=dgram_grouping,
                                dgram_function=dgram_function,
                                legend=legend,
                                xlim=xlim)

  cowplot::plot_grid(genes_plot,
                     trend_plot,
                     dgram_plot + ggplot2::theme(legend.position="bottom"),
                     ncol=1,
                     align="v")

}

plotDomainogram <- function(umi4c,
                            grouping="condition",
                            dgram_function="quotient", # or "difference"
                            legend=c("Treatment"="darkorange3",
                                     "Reference"="darkorchid3"),
                            xlim=NULL) {
  factor <- unique(colData[, grouping])

  ## Normalize dgrams
  dgram <- dgram(umi4c)

  for (i in 1:length(dgram)) {
    dgram[[i]] <- dgram[[i]]*assays(umi4c)$norm_mat[,i]
  }

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
                                    ymin=rev(scales), ymax=rev(scales)+1,
                                    fill=value)) +
    ggplot2::scale_fill_gradientn(colors=c(darken(legend[2], factor=10),
                                           legend[2], "white",
                                           legend[1],
                                           darken(legend[1], factor=10)),
                                  na.value = NA,
                                  name=as.expression(bquote(.(lab_legend)*log[2]*" UMIs")),
                                  breaks=scales::pretty_breaks(n=4),
                                  guide = ggplot2::guide_colorbar(direction = "horizontal",
                                                                  title.position="top",
                                                                  barwidth=8)) +
    ggplot2::coord_cartesian(xlim=xlim)

  return(dgram_plot)
}


plotTrend <- function(umi4c,
                      grouping=c("condition", "replicate"),
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

  trend_df$grouping_var <- do.call(paste, trend_df[,grouping])

  trend_plot <-
    ggplot2::ggplot(trend_df) +
    ggplot2::geom_ribbon(ggplot2::aes(geo_coord, ymin=trend-sd, ymax=trend+sd,
                                      group=grouping_var),
                         fill="grey", alpha=0.8, color=NA) +
    ggplot2::geom_line(ggplot2::aes(geo_coord, trend,
                                    group=grouping_var,
                                    color=grouping_var)) +
    ggplot2::coord_cartesian(xlim=xlim, ylim=ylim) +
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
                                 genes_uni[,c(8,12)])


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
