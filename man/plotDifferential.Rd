% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotUMI4C.R
\name{plotDifferential}
\alias{plotDifferential}
\title{Plot differential contacts}
\usage{
plotDifferential(umi4c, grouping = NULL, colors = NULL, xlim = NULL)
}
\arguments{
\item{umi4c}{\linkS4class{UMI4C} object as generated by \code{\link{makeUMI4C}}.}

\item{grouping}{Grouping used for the different samples. If none available or 
want to add new groupings, run \code{\link{addGrouping}}.}

\item{colors}{Named vector of colors to use for plotting for each group.}

\item{xlim}{Limits for the plot x axis (genomic coordinates).}
}
\value{
Produces a plot of the fold changes at the differential regions
analyzed ghat are contained in the \linkS4class{UMI4C} object.
}
\description{
Plot differential contacts
}
\examples{
data("ex_ciita_umi4c")
ex_ciita_umi4c <- addGrouping(ex_ciita_umi4c, grouping="condition")

enh <- GRanges(c("chr16:10925006-10928900", "chr16:11102721-11103700"))
umi_dif <- fisherUMI4C(ex_ciita_umi4c, query_regions = enh, 
                       filter_low = 20, resize = 5e3)
plotDifferential(umi_dif)
}
