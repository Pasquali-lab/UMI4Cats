% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotUMI4C.R
\name{plotUMI4C}
\alias{plotUMI4C}
\title{Plot UMI4C data}
\usage{
plotUMI4C(
  umi4c,
  grouping = "condition",
  dgram_function = "quotient",
  dgram_plot = TRUE,
  colors = NULL,
  xlim = NULL,
  ylim = NULL,
  TxDb = TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene,
  longest = TRUE,
  rel_heights = c(0.25, 0.4, 0.12, 0.23),
  font_size = 14
)
}
\arguments{
\item{umi4c}{\linkS4class{UMI4C} object as generated by \code{\link{makeUMI4C}}.}

\item{grouping}{Grouping used for the different samples. If none available or 
want to add new groupings, run \code{\link{addGrouping}}.}

\item{dgram_function}{Function used for calculating the fold-change in the
domainogram plot, either "difference" or "quotient". Default: "quotient".}

\item{dgram_plot}{Logical indicating whether to plot the domainogram. If the
\linkS4class{UMI4C} object only contains one sample will be automatically
set to FALSE. Default: TRUE.}

\item{colors}{Named vector of colors to use for plotting for each group.}

\item{xlim}{Limits for the plot x axis (genomic coordinates).}

\item{ylim}{Limits of the trend y axis.}

\item{TxDb}{TxDb object to use for drawing the genomic annotation.}

\item{longest}{Logical indicating whether to plot only the longest
transcripts for each gene in the gene annotation plot.}

\item{rel_heights}{Numeric vector of length 3 or 4 (if differential plot)
indicating the relative heights of each part of the UMI4C plot.}

\item{font_size}{Base font size to use for the UMI4C plot. Default: 14.}
}
\value{
Produces a summary plot with all the information contained in the
UMI4C opject.
}
\description{
Produce a UMI-4C data plot containing the genes in the region, the
adaptative smoothen trend and the domainogram.
}
\examples{
data("ex_ciita_umi4c")
ex_ciita_umi4c <- addGrouping(ex_ciita_umi4c, grouping="condition")

plotUMI4C(ex_ciita_umi4c,
    dgram_plot = FALSE
)
}
