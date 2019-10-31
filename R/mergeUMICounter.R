#' umi4Cats mergeUMICounter
#'
#'@description
#'Merge UMI4c counter with genome track positions.
#'
#'@usage
#'mergeUMICounter(genomic_track, wk_dir)
#'
#'@inheritParams umi4CatsContacts
#'
#'@examples
#'\dontrun{
#'wk_dir = '/imppc/labs/lplab/share/marc/epimutations/processed/prove/MLH1_ctrl_umi4cats_python'
#'genomic_track = '/imppc/labs/lplab/share/marc/epimutations/processed/genomicTracks/genomic_tracks_hg19/dpnII_genomicTrack'
#'
#'mergeUMICounter(genomic_track, wk_dir)
#'}
#'
#'@export

mergeUMICounter <- function(genomic_track,
                            wk_dir){

  bowtie2 <- 'bowtie2'

  counts_path <- file.path(wk_dir, 'rst')

  umi4cCountsFiles <- list.files(counts_path,
                                 full.names = T,
                                 pattern = '\\_umi_counts.tsv$')


  counts_path <- file.path(wk_dir, 'rst')

  for(umi4cCounts in umi4cCountsFiles){

    fileName <- gsub("\\..*$", "", basename(umi4cCounts))

    dfUmi4cCounts <- read.table(umi4cCounts)

    dfUmi4cCounts <- data.frame(lapply(dfUmi4cCounts, as.character), stringsAsFactors=FALSE)

    colnames(dfUmi4cCounts) <- c('chr_bait', 'pos_bait', 'chr_contact', 'pos_contact', 'UMIs')

    dfGenomicTrack <- read.table(genomic_track, stringsAsFactors = F)

    # get coordinates of viewpoint using bowtie2
    viewpoint <- paste0(bait_seq, bait_pad, res_e)

    index <- gsub('\\..*$', '', ref_gen)

    viewpointPos <- system(paste("bowtie2 --quiet",
                                 "-x", index,
                                 "-c", viewpoint,
                                 "-N 0",
                                 "| samtools view",
                                 "| awk \'{print $3,$4}\'"),
                           intern = T)

    viewpointPos <- unlist(strsplit(viewpointPos, " "))
    chrVP <- viewpointPos[1]
    startVP <- as.numeric(viewpointPos[2])

    start10M <- startVP - 10e6
    end10M <- startVP + 10e6

    subDfGenomicTrack <- dfGenomicTrack[(dfGenomicTrack$V1 == chrVP) &
                                          (dfGenomicTrack$V2 >= start10M) &
                                          (dfGenomicTrack$V2 <= end10M),]

    subDfGenomicTrack <- subDfGenomicTrack[c(1,2)]

    subDfGenomicTrack$chr_bait <- chrVP
    subDfGenomicTrack$pos_bait <- startVP

    colnames(subDfGenomicTrack)[1:2] <- c('chr_contact', 'pos_contact')

    subDfGenomicTrack <- data.frame(lapply(subDfGenomicTrack, as.character), stringsAsFactors=FALSE)

    dfUmi4cCounts10M <- dplyr::left_join(subDfGenomicTrack,dfUmi4cCounts)
    dfUmi4cCounts10M$UMIs[is.na(dfUmi4cCounts10M$UMIs)] <- 0

    dfUmi4cCounts10M <- dfUmi4cCounts10M[,c('chr_bait', 'pos_bait', 'chr_contact', 'pos_contact', 'UMIs')]

    umiOutput <- file.path(wk_dir, 'rst', paste0(fileName, '_counts10M.tsv'))

    write.table(x = dfUmi4cCounts10M,
                file = umiOutput,
                row.names = F,
                quote = F,
                sep = '\t')
  }
}

