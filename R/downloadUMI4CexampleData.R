#' Download UMI4Cats example datasets
#'
#' Downloads the required UMI4Cats example datasets.
#' @param out_dir Output directory for the datasets, defaults to tempdir().
#' @param verbose Whether to print verbose messages or not. Default: TRUE.
#' @param reduced Whether to use a reduced dataset to make test functions run faster.
#' @return It creates the \code{output_dir} with the example UMI-4C files used
#'  by the vignette and examples. Takes advantage of the BiocFileCache package to
#'  make sure that the file has not been previously downloaded by the user.
#' @examples
#' \dontrun{
#'  # Using reduced data data to make example faster.
#' # Remove reduced=TRUE or set to FALSE to
#' # download the full dataset.
#'
#' path <- downloadUMI4CexampleData(reduced = TRUE)
#' }
#' @import BiocFileCache
#' @importFrom utils download.file untar
#' @export
downloadUMI4CexampleData <- function(out_dir = tempdir(),
                                     verbose = TRUE,
                                     reduced = FALSE) {

    if (reduced) {
        file_url <- "http://gattaca.imppc.org/genome_browser/lplab/UMI4Cats_data_reduced.tar.gz" #https://ndownloader.figshare.com/files/25142468"
        rname <- "UMI4Cats_data_reduced"
    } else {
        file_url <- "http://gattaca.imppc.org/genome_browser/lplab/UMI4Cats_data.tar.gz" #"https://ndownloader.figshare.com/files/24309458"
        rname <- "UMI4Cats_data"
    }

    bfc <- .getCache()
    rid <- bfcquery(bfc, rname, "rname", exact = TRUE)$rid

    ## Add resource if it isn't cached
    if (!length(rid)) {
        if (verbose) message( "Downloading UMI4Cats data" )

        rid <- names(bfcadd(bfc, rname, file_url))
    }

    ## Update resource
    if (!isFALSE(bfcneedsupdate(bfc, rid))) {
        bfcdownload(bfc, rid)
    }

    untar(bfcrpath(bfc, rids = rid), verbose=TRUE, exdir=file.path(out_dir))

    return(file.path(out_dir, rname))
}


#' Get BiocFileCache object
#' @return Returns BFC object with the cache for the UMI4Cats package
.getCache <-
    function()
    {
        cache <- rappdirs::user_cache_dir(appname="UMI4Cats")
        BiocFileCache::BiocFileCache(cache)
    }
