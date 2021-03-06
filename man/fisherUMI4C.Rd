% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/differentialFisherUMI4C.R
\name{fisherUMI4C}
\alias{fisherUMI4C}
\title{Differential UMI4C contacts using Fisher's Exact test}
\usage{
fisherUMI4C(
  umi4c,
  grouping = "condition",
  query_regions,
  resize = NULL,
  window_size = 5000,
  filter_low = 50,
  padj_method = "fdr",
  padj_threshold = 0.05
)
}
\arguments{
\item{umi4c}{UMI4C object as generated by \code{makeUMI4C} or the
\code{UMI4C} constructor.}

\item{grouping}{Name of the column in colData used to merge the samples or
replicates. If none available or want to add new groupings, run 
\code{\link{addGrouping}}. Default: "condition".}

\item{query_regions}{\code{GenomicRanges} object or \code{data.frame}
containing the region coordinates used to perform the differential analysis.}

\item{resize}{Width in base pairs for resizing the \code{query_regions}.
Default: no resizing.}

\item{window_size}{If \code{query_regions} are not defined, wil bin region in
\code{window_size} bp and perform the analysis using this windows.}

\item{filter_low}{Either the minimum median UMIs requiered to perform
Fisher's Exact test or \code{FALSE} for performing the test in all windows.}

\item{padj_method}{Method for adjusting p-values. See
\code{\link[stats]{p.adjust}} for the different methods.}

\item{padj_threshold}{Numeric indicating the adjusted p-value threshold to
use to define significant differential contacts.}
}
\value{
Calculates statistical differences between UMI-4C experiments.
}
\description{
Using the UMIs inside \code{query_regions} performs Fisher's Exact test to
calculate significant differences between contact intensities.
}
\details{
This function calculates the
 overlap of fragment ends with either the provided \code{query_regions} or
 the binned region using \code{window_size}. The resulting number of UMIs in
 each \code{query_region} will be the \emph{sum} of UMIs in all overlapping
 fragment ends. As a default, will filter out those regions whose median
 UMIs are lower than \code{filter_low}.

Finally, a contingency table for each \code{query_reegions} or \code{window}
that passed the \code{filter_low} filter is created as follows:
\tabular{rcc}{
     \tab \emph{query_region} \tab \emph{region}\cr
    \emph{Reference} \tab n1 \tab N1-n1\cr
    \emph{Condition} \tab n2 \tab N2-n2
    }
and the Fisher's Exact test is performed. Obtained p-values are then
converted to adjusted p-values using \code{padj_method} and the results list
is added to the \code{UMI4C} object.
}
\examples{
data("ex_ciita_umi4c")
ex_ciita_umi4c <- addGrouping(ex_ciita_umi4c, grouping="condition")

# Perform differential test
enh <- GRanges(c("chr16:10925006-10928900", "chr16:11102721-11103700"))
umi_dif <- fisherUMI4C(ex_ciita_umi4c, query_regions = enh, 
                       filter_low = 20, resize = 5e3)
resultsUMI4C(umi_dif)                        
}
