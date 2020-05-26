#' UMI4Cats: A package for analyzing UMI-4C chromatin contact data
#'
#' The UMI4Cats package provides functions for the pre-processing, analysis and
#' visualization of UMI-4C chromatin contact data.
#'
#' @section File preparation:
#' There are two different functions that can be used to prepare the files
#' for analyzing them with UMI4Cats:
#' \enumerate{
#'     \item \code{\link{demultiplexFastq}}. Demultiplex reads belonging to
#'     different viewpoints from a paired-end FastQ file.
#'     \item \code{\link{digestGenome}}. Digest the reference genome of choice
#'     using a given restriction sequence.
#' }
#'
#' @section Processing:
#' The pre-processing functions are wrapped in the \code{\link{contactsUMI4C}}
#' main function. This function will sequentially run the following steps:
#' \enumerate{
#'     \item \code{\link{prepUMI4C}}. Filter specific and high quality reads.
#'     \item \code{\link{splitUMI4C}}. Split reads by the restriction sequence.
#'     \item \code{\link{alignmentUMI4C}}. Align reads to the reference genome.
#'     \item \code{\link{counterUMI4C}}. Apply UMI counting algorithm to quantify
#'     the interactions with the viewpoint.
#' }
#'
#' The statistics from the samples analyzed with the \code{\link{contactsUMI4C}}
#' function can be extracted and visualized with the function
#' \code{\link{statsUMI4C}}.
#'
#' @section Analysis:
#' The analysis of UMI-4C data is wrapped in the construction of an object of
#' \linkS4class{UMI4C} class by the creator function \code{\link{makeUMI4C}}.
#' This function will group your samples according to the variable you provided
#' in the \code{grouping} argument (default: "condition") and then normalize it
#' to \code{ref_umi4c}.
#'
#' The differential analysis is performed with the function \code{\link{fisherUMI4C}},
#' which will return a \linkS4class{UMI4C} object containing the results of the
#' differential test. You can access these results with the method \code{\link{results}}.
#'
#' @section Visualization:
#' An integrative plot showing the results stored inside the \linkS4class{UMI4C}
#' object can be generated with the function \code{\link{plotUMI4C}}.
#'
#' @docType package
#' @name UMI4Cats
NULL

utils::globalVariables(c(
    "factors", "scales", "value",
    "hg19_gene_annoation_ensemblv75", "stepping",
    "gene_name", "geo_coord", "grouping_var",
    "relative_position", "sample_id", "al_mapped",
    "al_unmapped", "variable",
    "queryHits", "subjectHits", "subsetByOverlaps",
    "assay", "metadata", "rowRanges", "assays", "assays<-",
    "assay", "colData", "IRanges", "SimpleList",
    "metadata<-", "rowRanges<-", "UMIs"
))
