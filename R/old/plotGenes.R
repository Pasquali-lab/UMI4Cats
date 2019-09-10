#' Plot genes
#'
#' @param window GRanges object with coordinates from which to plot existing genes.
#' @param protein_coding Plot only protein coding genes. Default: TRUE
#' @export
#' @importFrom magrittr %>%
plotGenes <- function(window,
                      protein_coding=TRUE) {
  xlim <- c(GenomicRanges::start(window),
            GenomicRanges::end(window))

  ## Get gene names in region
  genes.sel <- unique(IRanges::subsetByOverlaps(genes, window))
  if(protein_coding) genes.sel <- genes.sel[genes.sel$gene_biotype=="protein_coding",]

  ## Add stepping
  genes.step <- addStepping(genes.sel,
                            window,
                            2)
  genes.df <- data.frame(genes.step)

  ## Get exons
  exons <- getGenesAndExons(unique(genes.sel$external_gene_name))

  ## Edit genes
  genes.df <- dplyr::left_join(genes.df, exons)
  distance <- GenomicRanges::width(window)*0.01
  genes.df$modEnd <- genes.df$transcript_end + distance

  genes.df.uni <- unique(genes.df[,!grepl("exon", colnames(genes.df))])

  ## Edit starts and ends
  genes.df.uni$transcript_start[genes.df.uni$transcript_start < GenomicRanges::start(window)] <- GenomicRanges::start(window)
  genes.df.uni$transcript_end[genes.df.uni$transcript_end > GenomicRanges::end(window)] <- GenomicRanges::end(window)

  ## Plot genes--------------
  genesPlot <-
    ggplot2::ggplot(data=genes.df.uni) +
    ggplot2::geom_segment(data=genes.df.uni,
                          ggplot2::aes(x=transcript_start, y=stepping,
                                       xend=transcript_end, yend=stepping)) +
    ggplot2::geom_rect(data=genes.df,
                       ggplot2::aes(xmin=exon_chrom_start, xmax=exon_chrom_end,
                                    ymin=(stepping-0.3), ymax=(stepping+0.3)),
                                    fill="grey39", color="grey39") +
    ggplot2::geom_text(data=genes.df.uni,
                       ggplot2::aes(x=transcript_end, y=stepping),
                       label=as.character(genes.df.uni$external_gene_name),
                       colour="black",
                       hjust=0, fontface=3, nudge_x=distance,
                       size=3) +
    scaleXCoordinates(chr=as.character(GenomicRanges::seqnames(window)),
                      limits=xlim,
                      breaks=scales::pretty_breaks(n=4),
                      expand=c(0,0)) +
    themeYblank()

  return(genesPlot)

}

getGenesAndExons <- function(external_gene_name) {

  ## Get hg19 mart
  ensembl <-  biomaRt::useEnsembl(biomart="ensembl", GRCh=37,
                                dataset="hsapiens_gene_ensembl")

  ## Obtain all genes and transcripts
  exons <- biomaRt::getBM(attributes=c("exon_chrom_start", "exon_chrom_end", "ensembl_gene_id", "external_gene_name",
                                       "transcript_length", "transcript_start", "transcript_end",
                                       "ensembl_transcript_id"),
                          filters="external_gene_name",
                          values=external_gene_name,
                          mart=ensembl)

  longest <- exons %>%
    dplyr::group_by(external_gene_name) %>%
    dplyr::summarise(max=unique(ensembl_transcript_id[which(transcript_length==max(transcript_length))])[1])

  exons <- exons[exons$ensembl_transcript_id %in% longest$max,]

  return(exons)
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

