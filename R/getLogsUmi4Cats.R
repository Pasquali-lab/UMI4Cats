# getLogsUmi4Cats 
#'
#'@description
#' Plots stats of UMI4Cats filtering (based on the presence of bait and restriction enzyme) and alignment process. 
#' Results will be save in wk_dir/rst.
#'
#'@inheritParams umi4CatsContacts
#'
#'@examples
#'\dontrun{
#'raw_dir = '/imppc/labs/lplab/share/marc/epimutations/processed/prove/MLH1_ctrl_umi4cats_python'
#'wk_dir = '/imppc/labs/lplab/share/marc/epimutations/raw'
#'threads = '2'
#'
#'getLogsUmi4Cats(raw_dir = raw_dir, wk_dir = wk_dir, threads = threads)
#'}
#'
#'@export

getLogsUmi4Cats <- function(raw_dir,
                            wk_dir,
                            threads){
  
  # define variables
  samtools <- 'samtools'
  outPath <- file.path(wk_dir, 'rst')
  
  # parse files 
  raw_files <- list.files(raw_dir,
                          pattern = '_R1.fastq')
  
  raw_files <- gsub('_R1.fastq', '', raw_files)
  
  wk_files <- list.files(file.path(wk_dir, 'prep'),
                         pattern = '_prefiltered_R1.fastq')
  
  wk_files <- gsub('_prefiltered_R1.fastq', '', wk_files)
  
  samples <- unique(raw_files, wk_files)
  
  dfCounts <- data.frame()
  
  # get filtering and alignment stats and save as df
  for(sample in samples){
    
    # filter log
    raw_file <- list.files(raw_dir,
                           pattern = paste0(sample, '_R1.fastq'),
                           full.names = T)
    
    wk_file <- list.files(file.path(wk_dir, 'prep'),
                          pattern = paste0(sample, '_prefiltered_R1.fastq'),
                          full.names = T)
    
    sam_file <- list.files(file.path(wk_dir, 'alignment'),
                           pattern = paste0(sample, '\\.sam$'),
                           full.names = T)
    
    rawCounts <- system(paste('wc -l', raw_file),
                        intern = T)
    
    rawCounts <- unlist(strsplit(rawCounts, " "))[1]
    
    wkCounts <- system(paste('wc -l', wk_file),
                       intern = T)
    
    wkCounts <- unlist(strsplit(wkCounts, " "))[1]
    
    dfTmp <- t(data.frame(c(sample, rawCounts, wkCounts), stringsAsFactors = F))
    
    colnames(dfTmp) <- c('sample', 'raw', 'passFilter')
    
    # alignment log
    samFlagstat <- system(paste(samtools,
                                'flagstat',
                                "-@", threads,
                                sam_file),
                          intern = T)
    
    samFlagstat <- unlist(lapply(samFlagstat, function(x) unlist(strsplit(x, " "))[1]))
    
    samFlagstat <- t(data.frame(samFlagstat))
    
    samFlagstat <- data.frame(samFlagstat)
    
    colnames(samFlagstat) <- c('total', 'secondary', 'supplementary',
                               'duplicates', 'mapped', 'paired in sequencing',
                               'read1', 'read2', 'properly paired',
                               'with itself and mate mapped', 'singletons', 'with mate mapped to a different chr',
                               'with mate mapped to a different chr (mapQ>=5)')
    
    row.names(samFlagstat) <- seq(length(samFlagstat$total))
    
    dfTmp <- cbind(dfTmp, samFlagstat)
    
    dfCounts <- rbind(dfCounts, dfTmp)
    
  }

  rownames(dfCounts) <- seq(length(dfCounts$sample))
  
  # save logs
  write.table(dfCounts,
              file.path(outPath, "logs.txt"),
              row.names = F,
              quote = F,
              sep = "\t")
  
  # generate filter stats and plot
  dfCounts$raw <- as.numeric(levels(dfCounts$raw))
  dfCounts$passFilter <- as.numeric(levels(dfCounts$passFilter))
  
  dfCounts$nonPassFilter <- (dfCounts$raw - dfCounts$passFilter)/dfCounts$raw * 100
  
  dfCounts$passFilter <- dfCounts$passFilter/dfCounts$raw  * 100
  
  dfPass <- dfCounts[c('sample', 'passFilter', 'nonPassFilter')]
  
  dfPass <- reshape2::melt(dfPass, 'sample')
  
  dfPass <- dfPass[order(dfPass$sample, dfPass$value, decreasing = F),]
  
  gFiltered <- ggplot2::ggplot(dfPass, ggplot2::aes(x = sample, y = value, fill = variable)) + 
    ggplot2::geom_bar(stat = 'identity') +
    ggplot2::xlab('Sample') + ggplot2::ylab('% of reads') +
    ggplot2::geom_text(label = paste0(round(dfPass$value, 2),"%"), 
                       position=ggplot2::position_stack(0.5),) +
    ggplot2::scale_fill_discrete(name = "", labels = c("Filter passed ", "Filter not passed")) +
    cowplot::theme_cowplot() +
    ggplot2::theme(legend.position="top",
                   axis.text.x = ggplot2::element_text(angle = 90, hjust = 1))
  
  # generate alignment stats and plot
  dfSam <- dfCounts[c('sample', 'total', 'mapped')]
  
  dfSam$total <- as.numeric(levels(dfSam$total))
  dfSam$mapped <- as.numeric(levels(dfSam$mapped))
  
  dfSam$notMapped <- dfSam$total - dfSam$mapped
  dfSam$mapped <- (dfSam$mapped/dfSam$total)  * 100
  dfSam$notMapped <- dfSam$notMapped/dfSam$total  * 100
  
  dfSam <- reshape2::melt(dfSam, c('sample', 'total'))
  
  dfSam <- dfSam[order(dfSam$sample, dfSam$value, decreasing = F),]
  
  gAligned <- ggplot2::ggplot(dfSam, ggplot2::aes(x = sample, y = value, fill = variable)) + 
    ggplot2::geom_bar(stat = 'identity') +
    ggplot2::xlab('Sample') +
    ggplot2::geom_text(label = paste0(round(dfSam$value, 2),"%"), 
                       position=ggplot2::position_stack(0.5),) +
    ggplot2::scale_fill_discrete(name = "", labels = c("Mapped", "Not mapped")) +
    cowplot::theme_cowplot()  +
    ggplot2::theme(axis.text.y=ggplot2::element_blank(),
                   axis.title.y=ggplot2::element_blank(),
                   legend.position="top",
                   axis.text.x = ggplot2::element_text(angle = 90, hjust = 1))
  
  # grid plots
  gridPlot <- cowplot::plot_grid(gFiltered, gAligned, ncol = 2, align = 'hv')

  # save plots
  ggplot2::ggsave(file.path(outPath, "gridPlot.pdf"), gridPlot)
  ggplot2::ggsave(file.path(outPath, "gridPlot.png"), gridPlot)

}

