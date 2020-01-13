#' @title UMI4C class methods
#' @docType methods
#'
#' @description This page contains a summary of the different methods used to access the information contained
#' inside the \code{UMI4C-class} object. See the details section for more information on the different accessors.
#' @name UMI4C-methods
#'
#' @param object a \code{UMI4C-class} object.
#'
#' @details There are several accessors to easily retrive information from a \code{UMI4C-class} object:
#' \itemize{
#'   \item \code{dgram}: Returns a named list with the output domainograms for each sample.
#'   \item \code{bait}: Returns a \linkS4class{GRanges} object with the position of the bait.
#'   \item \code{trend}: Returns a data.frame in long format with the values of the adapative smoothen trend.
#'   \item \code{results}: Returns a \linkS4class{GRanges} or data.frame with the results of the differential analysis.
#' }
#' @examples
#' \dontrun{
#' # Create example UMI4C-class object with helper function
#' umi <- makeUMI4Cexample()
#'
#' # Access the different information inside the UMI4C object
#' dgram(umi)
#' bait(umi)
#' trend(umi)
#'
#' # Perform differential test
#' umi <- fisherUMI4C(umi)
#' results(umi)
#' }
#' @seealso UMI4C, UMI4C-class
NULL
