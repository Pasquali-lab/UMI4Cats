waldUMI4C <- function(umi4c,
                      query_regions=NULL,
                      subset="overlap", #sum
                      design=~condition,
                      normalized = TRUE,
                      padj_method = "fdr",
                      padj_threshold = 0.05) {
  
  ## Obtain only fragments in query_regions, if available
  if(!is.null(query_regions)) {
    if (subset=="overlap") umi4c <- subsetByOverlaps(umi4c, query_regions)   # Option 1: by overlaps
    else if (subset=="sum") umi4c <- combineUMI4C(umi4c, query_regions)   # Option 2: sum
  }
  
  ## Convert to dds
  dds <- UMI4C2dds(umi4c=umi4c,
                   design=design)
  
  ## Run DESeq2
  dds <- DESeq2::DESeq(dds, test="Wald")
  
  ## Convert back to UMI4C
  umi4c <- dds2UMI4C(umi4c=umi4c,
                     dds=dds,
                     design=~condition,
                     normalized=normalized,
                     padj_method=padj_method,
                     padj_threshold=padj_threshold)
  
  return(umi4c)
}


combineUMI4C <- function(umi4c,
                         query_regions) {
  matrix <- assay(umi4c)
  rowranges <- rowRanges(umi4c)
  
  hits <- findOverlaps(rowranges, query_regions)
  
  # Change id for mcol 4
  ids <- split(mcols(rowranges)[queryHits(hits),1], subjectHits(hits))
  
  mat_sp <- lapply(ids, function(x) matrix[x,])
  mat_sum <- lapply(mat_sp, function(x) if(is.null(dim(x))) x else colSums(x))
  mat_final <- do.call(rbind, mat_sum)
  
  umi4c_comb <- UMI4C(colData=colData(umi4c),
                 rowRanges=unique(query_regions[subjectHits(hits)]),
                 assays=SimpleList(umi = mat_final),
                 metadata=metadata(umi4c))
  
  return(umi4c_comb)
}