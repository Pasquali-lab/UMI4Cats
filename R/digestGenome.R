#' Digest reference genome
#'
#' Performs an \emph{in silico} digestion of a given reference genome using a
#' given restriction enzyme sequence.
#'
#' @param res_enz Character containing the restriction enzyme sequence.
#' @param cut_pos Numeric indicating the nucleotide position where restriction
#' enzyme cuts (zero-based) (for example, for DpnII is 0).
#' @param name_RE Restriction enzyme name.
#' @param ref_gen A BSgenome object of the reference genome.
#' @param sel_chr Character vector indicating which chromosomes to select for
#' the digestion. Default: chr1-22, chrX, chrY.
#' @param out_path Output path where to save the genomic track. The default is a
#'  directory named \code{digested_genome/} created in your working directory. The rda
#' objects are saved in folder named by the \code{ref_gene}_\code{name_RE} in the
#' \code{out_path} folder.
#' @return Creates a rda file for every chromosome defined in \code{sel_chr}.
#' @examples
#' library(BSgenome.Hsapiens.UCSC.hg19)
#' ref_gen <- BSgenome.Hsapiens.UCSC.hg19
#'
#' hg19_dpnii <- digestGenome(
#'     res_enz = "GATC",
#'     cut_pos = 0,
#'     name_RE = "dpnII",
#'     sel_chr = "chr16", # Only in chr16 to reduce example running time
#'     ref_gen = ref_gen,
#'     out_path = file.path(tempdir(), "digested_genome/")
#' )
#' @import BSgenome
#' @export
digestGenome <- function(res_enz,
    cut_pos,
    name_RE,
    ref_gen,
    sel_chr = paste0("chr", c(seq_len(22), "X", "Y")),
    out_path = "digested_genome/") {
    message(paste(
        "Generating digested genome using:\n",
        "> Restriction enzyme sequence:", res_enz, "\n",
        "> Restriction enzyme cut position:", cut_pos, "\n",
        "> Restriction enzyme name:", name_RE, "\n",
        "> Reference genome:", GenomeInfoDb::bsgenomeName(ref_gen), "\n",
        "> Output path:", out_path
    ))

    # TODO: Include RE database, given RE name use info on cutting sequence.
    # TODO: Deal with restriction enzymes that have Ns or positions with multiple
    # nucleotide options

    # Create directory if it doesn't exist
    dir.create(out_path, showWarnings = FALSE)
    out_path <- file.path(out_path, paste0(
        GenomeInfoDb::bsgenomeName(ref_gen),
        "_", name_RE
    ))
    dir.create(out_path, showWarnings = FALSE)


    # Select levels
    lvls <- GenomeInfoDb::seqnames(ref_gen)

    if (!is.null(sel_chr)) lvls <- lvls[lvls %in% sel_chr]

    if (length(lvls) == 0) stop("Levels ins 'sel_chr' are not present in your 'ref_gen' object. Try changing 'sel_chr' or setting it to 'NULL'")

    # Identify the recongnition sites for each chromosomal entry
    id_num <- 0
    for (chr in lvls) {
        sel_gen <- ref_gen[[chr]]

        # Find pattern in chr
        matches <- Biostrings::matchPattern(res_enz, sel_gen)
        IRanges::start(matches) <- IRanges::start(matches) - 1 + (cut_pos)
        IRanges::end(matches) <- IRanges::end(matches) - (nchar(res_enz) - cut_pos)

        gaps <- Biostrings::gaps(matches, start = 1, end = length(sel_gen))
        IRanges::end(gaps) <- IRanges::end(gaps) + 1

        digested_genome_gr <- GenomicRanges::GRanges(
            seqnames = chr,
            ranges = IRanges::IRanges(gaps)
        )

        # Add unique identifier for fragments
        digested_genome_gr$id <- paste0(
            "fragment_", name_RE, "_",
            (1 + id_num):(id_num + length(digested_genome_gr))
        )
        id_num <- id_num + length(digested_genome_gr)

        # Save digested genome
        out_track <- file.path(out_path, paste0(chr, ".rda"))
        save(digested_genome_gr, file = out_track)
    }

    message(paste("Finished genome digestion."))

    # Return path of the digested genome invisibly
    invisible(out_path)
}
