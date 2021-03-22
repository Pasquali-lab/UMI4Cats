#' Create gene annotation object
#' @inheritParams plotGenes
#' @return GRanges object with the gene annotation in the window.
#' @examples
#' library(TxDb.Hsapiens.UCSC.hg19.knownGene)
#' window <- GRanges("chr16:11298262-11400036")
#' gene_anno <- createGeneAnnotation(
#'     window = window,
#'     TxDb = TxDb.Hsapiens.UCSC.hg19.knownGene
#' )
#' @import GenomicRanges
#' @import TxDb.Hsapiens.UCSC.hg19.knownGene
#' @import org.Hs.eg.db
#' @export
createGeneAnnotation <- function(window,
    TxDb = TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene,
    longest = TRUE) {
    trans <- GenomicFeatures::transcriptsByOverlaps(TxDb,
        ranges = window,
        columns = c(
            "tx_id", "tx_name",
            "gene_id"
        )
    )
    trans <- trans[vapply(trans$gene_id, length, FUN.VALUE = integer(1)) > 0, ]

    if (length(trans) == 0) {
        return(trans)
    } else {
        trans$gene_id <- unlist(trans$gene_id)

        ## Select longest transcripts
        if (longest) {
            trans <- trans[order(trans$gene_id), ]
            trans_sp <- split(trans, trans$gene_id)
            sel <- unlist(lapply(trans_sp, function(x) width(x) == max(width(x))))
            trans <- trans[sel, ]
        }

        # If some have same width, select first example
        reps <- table(trans$gene_id) > 1
        reps <- names(reps)[reps]

        if (length(reps) > 0) {
            trans_uni <- c(
                trans[which(rowSums(vapply(reps, grepl, trans$gene_id, FUN.VALUE = logical(length(trans$gene_id)))) == 0), ],
                trans[vapply(reps, function(x) grep(x, trans$gene_id)[1], FUN.VALUE = integer(1))]
            )
        } else {
            trans_uni <- trans
        }

        trans_uni$type <- "GENE"

        ## Retrieve exons for selected transcripts
        exons <- GenomicFeatures::exons(TxDb,
                                        columns = c("tx_id"),
                                        filter = list("tx_id" = trans_uni$tx_id)
        )

        exons <- .unlistExons(exons)
        exons <- exons[exons$tx_id %in% trans_uni$tx_id]
        mcols(exons) <- dplyr::left_join(data.frame(mcols(exons)),
                                         data.frame(mcols(trans_uni)))

        exons$type <- "EXON"

        trans <- c(trans_uni, exons)
        trans$gene_id <- as.character(trans$gene_id)
        trans$tx_id <- as.character(trans$tx_id)
        trans$tx_name <- as.character(trans$tx_name)

        ## Get gene names
        sym <- unlist(annotate::lookUp(unique(trans$gene_id),
            data = "org.Hs.eg",
            what = "SYMBOL",
            load = TRUE
        ))

        df_name <- data.frame(
            gene_id = names(sym),
            gene_name = sym,
            stringsAsFactors = FALSE
        )

        mcols(trans) <- dplyr::left_join(
            data.frame(mcols(trans)),
            df_name,
            by = "gene_id"
        )

        return(trans)
    }
}

.unlistExons <- function(exons) {
    length <- vapply(exons$tx_id, length, FUN.VALUE=integer(1))
    reps <- unlist(mapply(rep, seq_len(length(exons)), each=length))
    exons_unl <- exons[reps]
    exons_unl$tx_id <- unlist(exons$tx_id)
    return(exons_unl)
}
