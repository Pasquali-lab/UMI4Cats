% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/differentialNbinomWaldTestUMI4C.R
\name{differentialNbinomWaldTestUMI4C}
\alias{differentialNbinomWaldTestUMI4C}
\title{Differential UMI4C contacts using DESeq2 Wald Test}
\usage{
differentialNbinomWaldTestUMI4C(
  umi4c,
  design = ~condition,
  normalized = TRUE,
  padj_method = "fdr",
  query_regions = NULL,
  padj_threshold = 0.05,
  penalty = 0.1,
  alpha = 20
)
}
\arguments{
\item{umi4c}{UMI4C object as generated by \code{makeUMI4C} or the
\code{UMI4C} constructor.}

\item{design}{A \code{formula} or \code{matrix}. The formula expresses how
the counts for each fragment end depend on the variables in \code{colData}.
See  \code{\link[DESeq2]{DESeqDataSet}}.}

\item{normalized}{Logical indicating if the function should return normalized
or raw UMI counts.  Default: TRUE.}

\item{padj_method}{The method to use for adjusting p-values, see
\code{\link[stats]{p.adjust}}.  Default: fdr.}

\item{query_regions}{\code{GRanges} object or \code{data.frame} containing
the coordinates of the genomic regions you want to use to perform the
analysis in specific genomic intervals. Default: NULL.}

\item{padj_threshold}{Numeric indicating the adjusted p-value threshold to
use to define significant differential contacts.  Default: 0.05.}

\item{penalty}{Amount of smoothing to be applied to the estimated functional 
parameter.  Default: 0.1.}

\item{alpha}{Approximate number of fragments desired for every basis function
of the B-spline basis. \code{floor((max(number of fragments)) / alpha)} is 
passed to \code{create.bspline.basis} as nbasis argument. 4 is the minimum
allowed value. Default: 20.}
}
\value{
\code{UMI4C} object with the DESeq2 Wald Test results.
}
\description{
Using a \code{UMI4C} object, infers the differences between conditions 
specified in \code{design} of the smooth monotone fitted values using a 
Wald Test from \code{DESeq2} package.
}
\details{
This function fits the signal trend of a variance stabilized count 
values using a symmetric monotone fit for the distance dependency. Then 
scales the raw counts across the samples to obtain normalized factors. 
Finally, it detects differences between conditions applying the DESeq2 Wald 
Test.
}
\examples{
\dontrun{
 files <- list.files(system.file("extdata", "CIITA", "count", package="UMI4Cats"),
                     pattern = "*_counts.tsv.gz",
                     full.names = TRUE
 )
# Create colData including all relevant information
colData <- data.frame(
  sampleID = gsub("_counts.tsv.gz", "", basename(files)),
  file = files,
  stringsAsFactors = FALSE
)

library(tidyr)
colData <- colData |>
  separate(sampleID,
           into = c("condition", "replicate", "viewpoint"),
           remove = FALSE
  )
  
# Make UMI-4C object including grouping by condition
umi <- makeUMI4C(
  colData = colData,
  viewpoint_name = "CIITA",
  grouping = NULL,
  bait_expansion = 2e6
)

umi_wald <- differentialNbinomWaldTestUMI4C(umi4c=umi,
                                            design=~condition,
                                            alpha = 100)
}
}
