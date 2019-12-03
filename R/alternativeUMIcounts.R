altCount <- function(digested_genome_gr,
                  pos_viewpoint,
                  filtered_bam_R1,
                  filtered_bam_R2,
                  res_enz,
                  count_dir,
                  filter_bp=10e6){
  # get coordinates of viewpoint
  viewpoint <- GenomicRanges::GRanges(seqnames=pos_viewpoint[1],
                                      ranges=IRanges::IRanges(start=as.numeric(pos_viewpoint[2]),
                                                              end=as.numeric(pos_viewpoint[3])))
  viewpoint_filter <- GenomicRanges::resize(viewpoint, width=filter_bp*2, fix="center")
  
  #TODO: Digested genome as tabix or rda to load it faster
  # load digest genome, filter and transform to granges
  digested_genome_gr <- IRanges::subsetByOverlaps(digested_genome_gr, viewpoint_filter)
  
  # read bam and transform to a granges
  bam_R1_gr <- GenomicAlignments::readGAlignments(filtered_bam_R1,
                                                  use.names=TRUE,
                                                  param=Rsamtools::ScanBamParam(which=viewpoint_filter))
  bam_R2_gr <- GenomicAlignments::readGAlignments(filtered_bam_R2,
                                                  use.names=TRUE,
                                                  param=Rsamtools::ScanBamParam(which=viewpoint_filter))
  
  mcols(bam_R1_gr)$header <- names(bam_R1_gr)
  mcols(bam_R1_gr)$umi <- sapply(strsplit(names(bam_R1_gr), ":"),
                                 function(x) x[1])
  mcols(bam_R1_gr)$readID <- sapply(strsplit(names(bam_R1_gr), ":"),
                                 function(x) paste0(x[2], "_", x[3]))
  
  mcols(bam_R2_gr)$header <- names(bam_R2_gr)
  mcols(bam_R2_gr)$umi <- sapply(strsplit(names(bam_R2_gr), ":"),
                                 function(x) x[1])
  mcols(bam_R2_gr)$readID <- sapply(strsplit(names(bam_R2_gr), ":"),
                                    function(x) paste0(x[2], "_", x[3]))
  
  
  gr1 <- GenomicAlignments::granges(bam_R1_gr, use.mcols=TRUE, use.names=FALSE)
  gr2 <- GenomicAlignments::granges(bam_R2_gr, use.mcols=TRUE, use.names=FALSE)
  
  gr <- c(gr1, gr2)
  
  gr_sp <- GenomicRanges::GRangesList(GenomicRanges::split(gr, gr$readID))
  
  n <- elementNROWS(gr_sp)
  
  ## Regions containing only the viewpoint
  uni <- unlist(gr_sp[n==1])
  uni <- regioneR::joinRegions(uni)
  
  ligations <- gr_sp[n>1]
  ligations_noview <- IRanges::subsetByOverlaps(unlist(ligations), uni, invert=TRUE)
  lig_noview_unl <- GenomicRanges::GRangesList(GenomicRanges::split(ligations_noview, 
                                                                ligations_noview$readID))
  
  # table(elementNROWS(lig_noview_unl))
  # table(elementNROWS(reduce(lig_noview_unl)))
  
  ## Select representatives from each ligation
  reps <- unlist(lig_noview_unl[end(lig_noview_unl) == max(end(lig_noview_unl))])
  
  dup_positions <- duplicated(start(reps)) & duplicated(end(reps))
  
  reps_pos <- reps[!dup_positions]
  
  ## Compare UMIs
  collapsed_umis <- c()
  umi_list <- unique(Biostrings::DNAStringSet(reps_pos$umi))
  
  while (length(umi_list) > 0) {
    
    compared_umi <- umi_list[[1]]
    collapsed_umis <- c(collapsed_umis, as.character(compared_umi))
    matches <- Biostrings::vcountPattern(compared_umi, umi_list, max.mismatch = 2)
    umi_list <- umi_list[!as.logical(matches)]
    
  }
  
  reps_pos_umis <- reps_pos[reps_pos$umi %in% collapsed_umis]
  
  # Select unique ligations
  final_ligations <- unlist(ligations[names(ligations) %in% unique(reps_pos_umis$readID)])
  
  # Compare each range with fragment id for each ligation
  hits <- findOverlaps(final_ligations, 
                       digested_genome_gr,
                       minoverlap=nchar(res_enz)+1)
  
  final <- final_ligations[queryHits(hits)]
  final$fragID <- as.character(mcols(digested_genome_gr)[subjectHits(hits),1])
  table <- as.data.frame(table(final$fragID))
  
  umis_df <- data.frame(digested_genome_gr)[,c(1:2,6)]
  umis_df <- dplyr::left_join(umis_df, table, by=c(V4="Var1"))
  
  viewpoint <- data.frame(subsetByOverlaps(digested_genome_gr, uni))[,c(1,2)]
  
  final_umis <- cbind(viewpoint[rep(1, nrow(umis_df)),],
                      umis_df[,-3])
  colnames(final_umis) <- c('chr_bait', 'pos_bait',
                            'chr_contact', 'pos_contact',
                            'UMIs')
  final_umis$UMIs[is.na(final_umis$UMIs)] <- 0
  
  file_name <- strsplit(basename(filtered_bam_R1), "_R1")[[1]][1]
  counts_file <- file.path(count_dir, paste0(file_name, '_altCounts.tsv'))
  
  write.table(x = final_umis,
              file = counts_file,
              row.names = F,
              quote = F,
              sep = '\t')

}