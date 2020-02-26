#' Create gene annotation object
#' @inheritParams plotGenes
#' @return GRanges object with the gene annotation in the window.
#' @import GenomicRanges
createGeneAnnotation <- function(window,
                                 TxDb=TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene,
                                 longest=TRUE) {
  trans <- GenomicFeatures::transcriptsByOverlaps(TxDb,
                                                  ranges=window,
                                                  columns=c("tx_id", "tx_name", "gene_id"))
  trans <- trans[sapply(trans$gene_id, length)>0,]
  trans$gene_id <- unlist(trans$gene_id)

  ## Select longest transcripts
  if (longest) {
    trans <- trans[order(trans$gene_id),]
    trans_sp <- split(trans, trans$gene_id)
    sel <- unlist(lapply(trans_sp, function(x) width(x) == max(width(x))))
    trans <- trans[sel,]
  }

  # If some have same width, select first example
  reps <- table(trans$gene_id) > 1
  reps <- names(reps)[reps]

  trans_uni <- c(trans[which(rowSums(sapply(reps, grepl, trans$gene_id))==0),],
                 trans[sapply(reps, function(x) grep(x, trans$gene_id)[1])])

  trans_uni$type <- "GENE"

  ## Retrieve exons for selected transcripts
  exons <- GenomicFeatures::exons(TxDb,
                                  columns=c("tx_id", "tx_name", "gene_id"),
                                  filter=list("tx_id"=trans_uni$tx_id))
  exons$tx_id <- sapply(exons$tx_id, function(x) x[x %in% trans_uni$tx_id])
  exons$tx_name <- sapply(exons$tx_name, function(x) x[x %in% trans_uni$tx_name])
  exons$gene_id <- sapply(exons$gene_id, function(x) x[x %in% trans_uni$gene_id])

  exons$type <- "EXON"

  trans <- c(trans_uni, exons)

  ## Get gene names
  sym <- unlist(annotate::lookUp(unique(trans$gene_id),
                                 data='org.Hs.eg',
                                 what="SYMBOL",
                                 load=TRUE))

  df_name <- data.frame(gene_id=names(sym),
                        gene_name=sym,
                        stringsAsFactors = FALSE)

  mcols(trans) <- dplyr::left_join(data.frame(mcols(trans)),
                                   df_name)

  return(trans)
}


