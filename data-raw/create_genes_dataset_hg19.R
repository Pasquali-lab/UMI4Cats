ensembl <-  biomaRt::useEnsembl(biomart="ensembl", GRCh=37,
                                dataset="hsapiens_gene_ensembl")

genes <- biomaRt::getBM(attributes=c("chromosome_name", "start_position", "end_position",
                                     "ensembl_gene_id", "external_gene_name","gene_biotype"),
                        mart=ensembl)

genes <- regioneR::toGRanges(genes)

ensembldb::seqlevelsStyle(genes) <- "UCSC"
genes <- genes[GenomicRanges::seqnames(genes) %in% paste0("chr", c(1:22, "X", "Y")),]

usethis::use_data(genes, overwrite=T)
usethis::use_data_raw()
