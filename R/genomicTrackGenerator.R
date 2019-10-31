#' Genomic track generator
#'
#' In silico digestion of a reference genome
#' @param cs5p 5' restriction enzyme cut sequence. In the GA|TC case "GA".
#' @param cs3p 3' restriction enzyme cut sequence. In the GA|TC case "TC"
#' @param nameRe Restriction enzyme name .
#' @param refgen A BSgenome object of the reference genome.
#' @param outPath Output path where to save the genomic track.
#' @examples
#' \dontrun{
#' output <- '/imppc/labs/lplab/share/marc/refgen/genomicTracksR'
#' the selected RE in this case is DpnII (|GATC), so the cs5p is "" and cs3p is GATC
#' cs5p <- ""
#' cs3p <- "GATC"
#' nameRe <- 'dpnII'
#' refgen <- BSgenome.Hsapiens.UCSC.hg19
#'
#' genomicTrackGenerator <- function(cp5p, cs3p, nameRe, refgen, output)
#'}
#'
#' @export
genomicTrackGenerator <- function(cs5p,
                                  cs3p,
                                  nameRe,
                                  refgen,
                                  outPath){

  genomeTrack = data.frame()

  re <- paste0(cs5p, cs3p)
  cp5p <- as.numeric(width(cs5p))
  cp3p <- as.numeric(width(cs3p))

  # Identify the recongnition sites for each chromosomal entry

  # Generate a dataframe with the length of MspI digested fragments

  message(paste('Generating genome track using:\n',
                '5\' restriction enzyme cut sequence:', cs5p, '\n',
                '3\' restriction enzyme cut sequence:', cs3p, '\n',
                'Restriction enzyme name:', nameRe, '\n',
                'Reference genome:', bsgenomeName(refgen), '\n',
                'Output path:', outPath))

  for (chr in seq_along(refgen)){
    m <- Biostrings::matchPattern("GATC", refgen[[chr]])
    starts <- start(Biostrings::gaps(m))
    ends <- end(Biostrings::gaps(m))
    temp_df <- data.frame(start = starts - cp3p, end = ends + cp5p, chr = seqnames(refgen)[chr]) # add re nucleotide length
    temp_df$start[1] <- temp_df$start[1] + cp3p # correct first position chr
    temp_df <- temp_df[c("chr","start","end")]
    genomeTrack <- rbind(genomeTrack,temp_df)
  }

  # Save genomicTrack

  dir.create(outPath, showWarnings = FALSE) # Create directory if it doesn't exist
  outTrack <- file.path(outPath, paste0(nameRe, '_genomicTrack.tsv'))

  write.table(genomeTrack,
              outTrack,
              col.names = F,
              row.names = F,
              quote = F,
              sep = "\t")

  message(paste('Finished genome track generation using:\n',
                '5\' restriction enzyme cut sequence:', cs5p, '\n',
                '3\' restriction enzyme cut sequence:', cs3p, '\n',
                'Restriction enzyme name:', nameRe, '\n',
                'Reference genome:', bsgenomeName(refgen), '\n',
                'Output path:', outPath))
}
