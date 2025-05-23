#' UMI4C Contacts Processing
#'
#' Using demultiplexed FastQ files as input, performs all necessary steps to end
#'  up with a tsv file summarizing the
#' restriction enzyme fragments and the number of UMIs supporting that specific
#' contact with the viewpoint (bait) of interest.
#'
#' @param fastq_dir Path of the directory containing the FastQ files (compressed
#'  or uncompressed).
#' @param wk_dir Working directory where to save the outputs generated by the
#' UMI-4c analysis.
#' @param file_pattern Character that can be used to filter the files you want
#' to analyze in the \code{fastq_dir}.
#' @param bait_seq Character containing the bait primer sequence.
#' @param bait_pad Character containing the pad sequence (sequence between the
#' bait primer and the restriction enzyme sequence).
#' @param res_enz Character containing the restriction enzyme sequence.
#' @param cut_pos Numeric indicating the nucleotide position where restriction
#' enzyme cuts (zero-based) (for example, for DpnII is 0).
#' @param digested_genome Path for the digested genome file generated using the
#' \code{\link{digestGenome}} function.
#' @param ref_gen A BSgenome object of the reference genome.
#' @param sel_seqname A character with the chromosome name to focus the
#' search for the viewpoint sequence.
#' @param threads Number of threads to use in the analysis. Default=1.
#' @param numb_reads Number of lines from the FastQ file to load in each loop.
#' If having memory size problems, change it to a smaller number. Default=1e9.
#' @param rm_tmp Logical indicating whether to remove temporary files (sam and
#' intermediate bams). TRUE or FALSE. Default=TRUE.
#' @param min_flen Minimal fragment length to use for selecting the fragments.
#' Default=20
#' @param filter_bp Integer indicating the bp upstream and downstream of the
#' viewpoint to select for further analysis. Default=10e6
#' @param bowtie_index Path and prefix of the bowtie index to use for the
#' alignment.
#' @return This function is a combination of calls to other functions that
#' perform the necessary steps for processing
#' UMI-4C data.
#' @examples
#' if (interactive()) {
#' path <- downloadUMI4CexampleData()
#'
#' hg19_dpnii <- digestGenome(
#'     cut_pos = 0,
#'     res_enz = "GATC",
#'     name_RE = "DpnII",
#'     ref_gen = BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19,
#'     out_path = file.path(path, "digested_genome")
#' )
#'
#'
#' raw_dir <- file.path(path, "CIITA", "fastq")
#'
#' contactsUMI4C(
#'     fastq_dir = raw_dir,
#'     wk_dir = file.path(path, "CIITA"),
#'     bait_seq = "GGACAAGCTCCCTGCAACTCA",
#'     bait_pad = "GGACTTGCA",
#'     res_enz = "GATC",
#'     cut_pos = 0,
#'     digested_genome = hg19_dpnii,
#'     bowtie_index = file.path(path, "ref_genome", "ucsc.hg19.chr16"),
#'     threads = 1,
#'     numb_reads = 1e9,
#'     ref_gen = BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19,
#'     sel_seqname = "chr16"
#' )
#'
#' unlink(path, recursive=TRUE)
#' }
#' @export
contactsUMI4C <- function(fastq_dir,
    wk_dir,
    file_pattern = NULL,
    bait_seq,
    bait_pad,
    res_enz,
    cut_pos,
    digested_genome,
    bowtie_index,
    threads = 1,
    numb_reads = 1e9,
    rm_tmp = TRUE,
    min_flen = 20,
    filter_bp = 10e6,
    ref_gen,
    sel_seqname = NULL) {
    dir.create(wk_dir, showWarnings = FALSE) # Create working dir
    # cut_pos <- as.character(cut_pos) # convert to character

    # get coordinates of viewpoint using biostrings
    pos_viewpoint <- getViewpointCoordinates(
        bait_seq = bait_seq,
        bait_pad = bait_pad,
        res_enz = res_enz,
        ref_gen = ref_gen,
        sel_seqname = sel_seqname
    )
    # Run steps of UMI4C analysis
    prepUMI4C(
        fastq_dir = fastq_dir,
        wk_dir = wk_dir,
        file_pattern = file_pattern,
        bait_seq = bait_seq,
        bait_pad = bait_pad,
        res_enz = res_enz,
        numb_reads = numb_reads
    )

    splitUMI4C(
        wk_dir = wk_dir,
        res_enz = res_enz,
        cut_pos = cut_pos,
        numb_reads = numb_reads,
        min_flen = min_flen
    )

    alignmentUMI4C(
        wk_dir = wk_dir,
        pos_viewpoint = pos_viewpoint,
        bowtie_index = bowtie_index,
        threads = threads,
        filter_bp = filter_bp
    )

    counterUMI4C(
        wk_dir = wk_dir,
        pos_viewpoint = pos_viewpoint,
        res_enz = res_enz,
        digested_genome = digested_genome,
        filter_bp = filter_bp
    )

    # Remove unnecessary folders
    if (rm_tmp) {
        unlink(paste0(wk_dir, "/", c("align", "prep", "split")),
            recursive = TRUE
        )
    }
}

#' Prepare UMI4C data
#'
#' Prepare the FastQ files for the further analysis by selecting reads with bait and
#' adding the respective UMI identifier for each read in its header.
#' @inheritParams contactsUMI4C
#' @return Creates a compressed FASTQ file in \code{wk_dir/prep} named
#' \code{basename(fastq)).fq.gz}, containing the filtered reads with the UMI
#' sequence in the header. A log file with the statistics is also generated
#' in \code{wk_dir/logs} named \code{umi4c_stats.txt}.
#' @examples
#' if (interactive()) {
#' path <- downloadUMI4CexampleData(reduced = TRUE)
#' raw_dir <- file.path(path, "CIITA", "fastq")
#'
#' prepUMI4C(
#'     fastq_dir = raw_dir,
#'     wk_dir = file.path(path, "CIITA"),
#'     bait_seq = "GGACAAGCTCCCTGCAACTCA",
#'     bait_pad = "GGACTTGCA",
#'     res_enz = "GATC"
#' )
#' }
#' @seealso \code{\link{contactsUMI4C}}.
#' @export
prepUMI4C <- function(fastq_dir,
    wk_dir,
    file_pattern = NULL,
    bait_seq,
    bait_pad,
    res_enz,
    numb_reads = 1e9) {
    message(paste(
        paste0("\n[", Sys.time(), "]"),
        "Starting prepUMI4C using:\n",
        "> Fastq directory:\n", fastq_dir, "\n",
        "> Work directory:", wk_dir, "\n",
        "> Bait sequence:", bait_seq, "\n",
        "> Bait pad:", bait_pad, "\n",
        "> Restriction enzyme:", res_enz, "\n",
        "> Number of reads loaded into memory:", numb_reads
    ))

    # create directory
    prep_dir <- file.path(wk_dir, "prep")
    dir.create(prep_dir, showWarnings = FALSE)

    # define variables
    fastq_files <- list.files(fastq_dir,
        pattern = "\\.fastq$|\\.fq$|\\.fq.gz$|\\.fastq.gz$",
        full.names = TRUE
    )

    if (!is.null(file_pattern)) {
        fastq_files <- fastq_files[grep(
            file_pattern,
            fastq_files
        )]
    }

    if (length(fastq_files) < 2) stop(paste("Non paired-end FASTQ files with the extension _RX.fastq, _RX.fq, _RX.fq.gz
                                        or _RX.fastq.gz in"), fastq_dir)

    fastqR1_files <- fastq_files[grep("_R1", fastq_files)]
    fastqR2_files <- fastq_files[grep("_R2", fastq_files)]


    # apply main function to files
    stats <- lapply(
        seq_len(length(fastqR1_files)),
        function(i) {
            .singlePrepUMI4C(
                fq_R1 = fastqR1_files[i],
                fq_R2 = fastqR2_files[i],
                bait_seq = bait_seq,
                bait_pad = bait_pad,
                res_enz = res_enz,
                prep_dir = prep_dir,
                numb_reads = numb_reads
            )
        }
    )

    # create stats file and save
    stats <- do.call(rbind, stats)
    dir.create(file.path(wk_dir, "logs"), showWarnings = FALSE) # Create logs dir
    utils::write.table(stats,
        file = file.path(wk_dir, "logs", "umi4c_stats.txt"),
        row.names = FALSE,
        sep = "\t",
        quote = FALSE
    )
}

#' Prepar fastq files at a given barcode.
#'
#' @param fq_R1 Fastq file R1.
#' @param fq_R2 Fastq file R2.
#' @param prep_dir Prep directory.
#' @inheritParams contactsUMI4C
#' @return Creates a compressed FASTQ file in \code{wk_dir/prep} named
#' \code{basename(Fastq)).fq.gz}, containing the filtered reads with the UMI
#' sequence in the header. A data.frame object with the statisitics is also
#' returned.
.singlePrepUMI4C <- function(fq_R1,
    fq_R2,
    bait_seq,
    bait_pad,
    res_enz,
    prep_dir,
    numb_reads = 1e9) {
    stream1 <- ShortRead::FastqStreamer(fq_R1, n = numb_reads)
    stream2 <- ShortRead::FastqStreamer(fq_R2, n = numb_reads)

    # Check output fastq files
    prep_fastqR1 <- paste0(gsub("\\..*", "", basename(fq_R1)), ".fq.gz")
    prep_fastqR2 <- paste0(gsub("\\..*", "", basename(fq_R2)), ".fq.gz")

    if (file.exists(file.path(prep_dir, prep_fastqR1))) unlink(file.path(prep_dir, prep_fastqR1))
    if (file.exists(file.path(prep_dir, prep_fastqR2))) unlink(file.path(prep_dir, prep_fastqR2))

    # Initialize global variables
    total_reads <- 0
    specific_reads <- 0
    filtered_reads <- 0

    # filter reads that not present bait seq + bait pad + re
    barcode <- paste0(bait_seq, bait_pad, res_enz)

    repeat {
        reads_fqR1 <- ShortRead::yield(stream1)
        reads_fqR2 <- ShortRead::yield(stream2)

        if (length(reads_fqR1) == 0) break
        if (length(reads_fqR1) != length(reads_fqR2)) stop("Different number of reads in R1 vs R2")

        total_reads <- total_reads + length(reads_fqR1) # Save total reads

        # for cases when the bait is to far from restriction enzyme
        if (nchar(barcode) > unique(width(reads_fqR1))) {
            barcode <- substr(barcode, 1, unique(width(reads_fqR1)))
        }

        barcode_reads_fqR1 <- reads_fqR1[grepl(barcode, ShortRead::sread(reads_fqR1))]
        barcode_reads_fqR2 <- reads_fqR2[grepl(barcode, ShortRead::sread(reads_fqR1))]

        specific_reads <- specific_reads + length(barcode_reads_fqR1) # Save specific reads

        # filter reads with less than 20 phred score --------------------
        filter20phred <- lapply(as(
            Biostrings::PhredQuality(Biostrings::quality(barcode_reads_fqR1)),
            "IntegerList"
        ), mean) >= 20 &
            lapply(as(Biostrings::PhredQuality(Biostrings::quality(barcode_reads_fqR2)), "IntegerList"), mean) >= 20

        filtered_reads_fqR1 <- ShortRead::ShortReadQ(
            ShortRead::sread(barcode_reads_fqR1)[filter20phred],
            Biostrings::quality(barcode_reads_fqR1)[filter20phred],
            ShortRead::id(barcode_reads_fqR1)[filter20phred]
        )

        filtered_reads_fqR2 <- ShortRead::ShortReadQ(
            ShortRead::sread(barcode_reads_fqR2)[filter20phred],
            Biostrings::quality(barcode_reads_fqR2)[filter20phred],
            ShortRead::id(barcode_reads_fqR2)[filter20phred]
        )

        # insert umi identifier (10 first bp of R2) to header of both R1 R2 files -----------
        umis <- stringr::str_sub(as.character(ShortRead::sread(filtered_reads_fqR2)), start = 1, end = 10)

        new_id_R1 <- paste0(umis, ":", "UMI4C:", seq(filtered_reads + 1, filtered_reads + length(filtered_reads_fqR1)), ":R1")
        new_id_R2 <- paste0(umis, ":", "UMI4C:", seq(filtered_reads + 1, filtered_reads + length(filtered_reads_fqR2)), ":R2")
        
        filtered_reads <- filtered_reads + length(filtered_reads_fqR1) # Return num filtered reads

        umi_reads_fqR1 <- ShortRead::ShortReadQ(
            ShortRead::sread(filtered_reads_fqR1),
            Biostrings::quality(filtered_reads_fqR1),
            Biostrings::BStringSet(new_id_R1)
        )

        umi_reads_fqR2 <- ShortRead::ShortReadQ(
            ShortRead::sread(filtered_reads_fqR2),
            Biostrings::quality(filtered_reads_fqR2),
            Biostrings::BStringSet(new_id_R2)
        )

        # write output fastq files
        ShortRead::writeFastq(umi_reads_fqR1,
            file.path(prep_dir, prep_fastqR1),
            mode = "a"
        )
        ShortRead::writeFastq(umi_reads_fqR2,
            file.path(prep_dir, prep_fastqR2),
            mode = "a"
        )
    }

    message(
        paste0("[", Sys.time(), "] "),
        "Finished sample ", strsplit(basename(fq_R1), "_R.")[[1]][1]
    )

    on.exit(close(stream1))
    on.exit(close(stream2), add = TRUE)

    # Construct stats data.frame
    stats <- data.frame(
        sample_id = strsplit(basename(fq_R1), "_R1")[[1]][1],
        total_reads = total_reads,
        specific_reads = specific_reads,
        filtered_reads = filtered_reads,
        stringsAsFactors = FALSE
    )
    return(stats)
}

#' Split UMI4C reads
#'
#' Split the prepared reads using the restrition enzyme information.
#' @inheritParams contactsUMI4C
#' @return Creates a compressed FASTQ file in \code{wk_dir/split} named
#' \code{basename(fastq)).fq.gz}, containing the
#' split reads based on the restriction enzyme used.
#' @examples
#' if (interactive()) {
#' path <- downloadUMI4CexampleData(reduced = TRUE)
#'
#' splitUMI4C(
#'     wk_dir = file.path(path, "CIITA"),
#'     res_enz = "GATC",
#'     cut_pos = 0
#' )
#' }
#' @export
splitUMI4C <- function(wk_dir,
    res_enz,
    cut_pos,
    numb_reads = 1e9,
    min_flen = 20) {
    message(paste(
        paste0("\n[", Sys.time(), "]"),
        "Starting splitUMI4C using:\n",
        "> Work directory:", wk_dir, "\n",
        "> Cut position:", cut_pos, "\n",
        "> Restriction enzyme:", res_enz, "\n",
        "> Number of reads loaded into memory:", numb_reads
    ))

    # create directory
    prep_dir <- file.path(wk_dir, "prep")
    split_dir <- file.path(wk_dir, "split")

    dir.create(split_dir, showWarnings = FALSE)

    prep_files <- list.files(prep_dir,
        pattern = ".gz$",
        full.names = TRUE
    )

    if (length(prep_files) < 2) stop(paste("Non paired-end prep FASTQ files with the extension _RX.fastq.gz
                                       or _RX.fq.gz in"), prep_dir)

    prep_files_R1 <- prep_files[grep("_R1", prep_files)]
    prep_files_R2 <- prep_files[grep("_R2", prep_files)]

    # run main function
    lapply(prep_files_R1, .singleSplitUMI4C,
        res_enz = res_enz, cut_pos = cut_pos, split_dir = split_dir,
        numb_reads = numb_reads
    )

    # run main function
    lapply(prep_files_R2, .singleSplitUMI4C,
        res_enz = res_enz, cut_pos = (nchar(res_enz) - cut_pos), split_dir = split_dir,
        numb_reads = numb_reads
    )
}

#' Split fastq files at a given restriction site.
#' @param fastq_file Fastq file path.
#' @param split_dir Directory where to save split files.
#' @inheritParams contactsUMI4C
#' @return Creates a compressed FASTQ file in \code{wk_dir/split} named
#' \code{basename(fastq)).fq.gz}, containing the split reads based on the
#' restriction enzyme used.
.singleSplitUMI4C <- function(fastq_file,
    res_enz,
    cut_pos,
    split_dir,
    min_flen = 20,
    numb_reads = 1e9) {

    # Use stream
    stream <- ShortRead::FastqStreamer(fastq_file, n = numb_reads)
    on.exit(close(stream))

    # Remove file if already exists
    split_fastq_name <- paste0(
        gsub("\\..*$", "", basename(fastq_file)),
        ".fq.gz"
    )
    filename <- file.path(split_dir, split_fastq_name)

    # Remove file if it already exists to avoid appending new reads
    if (file.exists(filename)) unlink(filename)

    repeat {
        # define variables and create objects
        prep_reads <- ShortRead::yield(stream)
        if (length(prep_reads) == 0) break

        prep_dna_string <- ShortRead::sread(prep_reads)

        ids <- ShortRead::id(prep_reads)

        # Find matches for the re sequence
        matches <- Biostrings::vmatchPattern(res_enz, prep_dna_string)
        matches <- as(matches, "CompressedIRangesList")
        IRanges::start(matches) <- as(
            IRanges::start(matches) - 1,
            "CompressedIntegerList"
        )
        IRanges::end(matches) <- as(
            IRanges::end(matches) - (nchar(res_enz) - cut_pos),
            "CompressedIntegerList"
        )

        # workaround for obtaining the cut position
        gaps <- IRanges::gaps(matches,
            start = 1,
            end = unique(nchar(as.character(prep_dna_string)))
        )

        ids_sel <- gaps
        IRanges::start(ids_sel) <- as(1, "IntegerList")
        IRanges::end(ids_sel) <- as(nchar(as.character(ids)), "IntegerList")

        list_seqs <- Biostrings::extractAt(prep_dna_string, gaps)
        list_quals <- Biostrings::extractAt(Biostrings::quality(Biostrings::quality(prep_reads)), gaps)
        list_ids <- Biostrings::extractAt(ids, ids_sel)

        fastq_entry <- ShortRead::ShortReadQ()
        fastq_entry@sread <- unlist(list_seqs)
        fastq_entry@quality <- ShortRead::FastqQuality(unlist(list_quals))
        fastq_entry@id <- unlist(list_ids)

        # Remove reads shorter than minimum fragment length
        fastq_entry <- fastq_entry[ShortRead::width(fastq_entry) >= min_flen]

        ShortRead::writeFastq(fastq_entry,
            file = filename,
            mode = "a"
        )
    }

    message(
        paste0("[", Sys.time(), "] "),
        "Finished sample ", basename(fastq_file)
    )
}

#' UMI4C alignment
#'
#' Align split UMI-4C reads to a reference genome using Bowtie2.
#' @inheritParams contactsUMI4C
#' @param pos_viewpoint GRanges object containing the genomic position of the
#' viewpoint. It can be generated by \code{getViewpointCoordinates} function.
#' @return Creates a BAM file in \code{wk_dir/align} named
#' "\code{basename(fastq))_filtered.bam}", containing the
#' aligned filtered reads. The alignment log is also generated in
#' \code{wk_dir/logs} named "\code{umi4c_alignment_stats.txt}".
#' @examples
#' if (interactive()){
#' path <- downloadUMI4CexampleData(reduced = TRUE)
#' alignmentUMI4C(
#'     wk_dir = file.path(path, "CIITA"),
#'     pos_viewpoint = GenomicRanges::GRanges("chr16:10972515-10972548"),
#'     bowtie_index = file.path(path, "ref_genome", "ucsc.hg19.chr16")
#' )
#' }
#' @export
alignmentUMI4C <- function(wk_dir,
    pos_viewpoint,
    bowtie_index,
    threads = 1,
    filter_bp = 10e6) {
    message(paste(
        paste0("\n[", Sys.time(), "]"),
        "Starting alignmentUMI4C using:\n",
        "> Work directory:", wk_dir, "\n",
        "> Viewpoint position:", paste0(
            as.character(GenomeInfoDb::seqnames(pos_viewpoint)),
            ":",
            as.character(IRanges::ranges(pos_viewpoint))
        ), "\n",
        "> Reference genome:", bowtie_index, "\n",
        "> Number of threads:", threads
    ))


    if (length(pos_viewpoint) == 0) stop(paste("Define viewpoint position"))

    # align split files
    split_dir <- file.path(wk_dir, "split")
    align_dir <- file.path(wk_dir, "align")

    dir.create(align_dir, showWarnings = FALSE)


    gz_files <- list.files(split_dir,
        pattern = ".gz$",
        full.names = TRUE
    )

    if (length(gz_files) != 0) {
        lapply(gz_files, R.utils::gunzip, overwrite = TRUE)
    }

    split_files <- list.files(split_dir,
        pattern = "\\.fastq$|\\.fq$",
        full.names = TRUE
    )
    if (length(split_files) < 2) {
        stop(
            paste("No paired-end split FASTQ files with the extension _RX.fastq.gz
                                          ,_RX.fq.gz, _RX.fastq or _RX.fq in"),
            split_dir
        )
    }

    stats <- lapply(split_files,
        .singleAlignmentUMI4C,
        align_dir = align_dir,
        threads = threads,
        bowtie_index = bowtie_index,
        pos_viewpoint = pos_viewpoint
    )

    stats <- do.call(rbind, stats)
    utils::write.table(stats,
        file = file.path(wk_dir, "logs", "umi4c_alignment_stats.txt"),
        row.names = FALSE, sep = "\t", quote = FALSE
    )
}

#' Align split fastq file
#' @inheritParams contactsUMI4C
#' @param split_file Split fastq file to align.
#' @param align_dir Directory where to save aligned files.
#' @param pos_viewpoint GRanges object containing the genomic position of the
#' viewpoint.
#' @return Creates a BAM file in \code{wk_dir/align} named
#' "\code{basename(fastq))_filtered.bam}", containing the aligned filtered
#' reads. A data.frame object with the statisitics is also returned.
#' @inheritParams contactsUMI4C
.singleAlignmentUMI4C <- function(split_file,
    align_dir,
    threads = 1,
    bowtie_index,
    pos_viewpoint,
    filter_bp = 10e6) {
    split_name <- gsub("\\..*$", "", basename(split_file))
    sam <- file.path(align_dir, paste0(split_name, ".sam"))
    bam <- file.path(align_dir, paste0(split_name, ".bam"))
    filtered_tmp_bam <- file.path(align_dir, paste0(
        split_name,
        "_filtered_tmp.bam"
    ))
    filtered_bam <- file.path(align_dir, paste0(split_name, "_filtered.bam"))
    # TODO:create bowtie2 index if it does not exist
    # align using bowtie2
    suppressMessages(Rbowtie2::bowtie2(
        seq1 = split_file,
        bt2Index = bowtie_index,
        "--threads", threads,
        samOutput = sam,
        overwrite = TRUE
    ))

    # sam to bam
    Rsamtools::asBam(sam, overwrite = TRUE)

    # keep reads in a 10M window from viewpoint
    pos_filter <- GenomicRanges::resize(pos_viewpoint,
        fix = "center",
        width = filter_bp * 2
    )

    param_10M <- Rsamtools::ScanBamParam(which = pos_filter)
    Rsamtools::filterBam(bam, filtered_tmp_bam, param = param_10M)

    # filter reads with 42mapq at least
    filter_mapq <- S4Vectors::FilterRules(list(mapq_filter = function(x) x$mapq >= 30))
    Rsamtools::filterBam(filtered_tmp_bam,
        filtered_bam,
        filter = filter_mapq,
        param = Rsamtools::ScanBamParam(what = "mapq")
    )

    message(
        paste0("[", Sys.time(), "] "),
        "Finished sample ", basename(split_file)
    )

    # Obtain stats for alignment
    stats <- data.frame(
        sample_id = gsub(".sam", "", basename(sam)),
        al_mapped = .getSummaryBam(gsub(".sam", ".bam", sam),
            mapped = TRUE
        ),
        al_unmapped = .getSummaryBam(gsub(".sam", ".bam", sam),
            mapped = FALSE
        ),
        al_secondary = .getSummaryBam(gsub(".sam", ".bam", sam),
            mapped = TRUE, secondary = TRUE
        )
    )
    return(stats)
}


#' UMI counting
#'
#' Algorithm for counting and collapsing the number of UMIs supporting a
#' specific ligation.
#' @inheritParams contactsUMI4C
#' @param pos_viewpoint GRanges object containing the genomic position of the
#' viewpoint.
#' @param filter_bp Integer indicating the bp upstream and downstream of the
#' viewpoint to select for further analysis. Default=10e6.
#' @return Creates a compressed tab-delimited file in \code{wk_dir/count} named
#' "\code{basename(fastq) _counts.tsv.gz}", containing the
#' coordinates for the viewpoint fragment, contact fragment and the number of
#' UMIs detected in the ligation.
#' @examples
#' if (interactive()) {
#' path <- downloadUMI4CexampleData(reduced = TRUE)
#'
#' hg19_dpnii <- digestGenome(
#'     cut_pos = 0,
#'     res_enz = "GATC",
#'     name_RE = "DpnII",
#'     sel_chr = "chr16", # digest only chr16 to make example faster
#'     ref_gen = BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19,
#'     out_path = file.path(path, "digested_genome")
#' )
#'
#' viewpoint <- GenomicRanges::GRanges("chr16:10972515-10972548")
#'
#' counterUMI4C(
#'     wk_dir = file.path(path, "CIITA"),
#'     pos_viewpoint = viewpoint,
#'     res_enz = "GATC",
#'     digested_genome = hg19_dpnii
#' )
#'
#' }
#' @details For collapsing different molecules into the same UMI, takes into
#' account the ligation position and the number of UMI sequence mismatches.
#' @export
counterUMI4C <- function(wk_dir,
    pos_viewpoint,
    res_enz,
    digested_genome,
    filter_bp = 10e6) {
    message(paste(
        paste0("\n[", Sys.time(), "]"),
        "Starting counterUMI4C using:\n",
        "> Work directory:", wk_dir, "\n",
        "> Viewpoint position:", paste0(
            as.character(GenomeInfoDb::seqnames(pos_viewpoint)),
            ":",
            as.character(IRanges::ranges(pos_viewpoint))
        ), "\n",
        "> Restriction enzyme:", res_enz, "\n",
        "> Digested genome:", digested_genome
    ))

    if (!exists("pos_viewpoint")) stop(paste("Define viewpoint position"))

    align_dir <- file.path(wk_dir, "align")
    count_dir <- file.path(wk_dir, "count")
    dir.create(count_dir, showWarnings = FALSE)

    # define variables
    aligned_files <- list.files(align_dir,
        pattern = "_filtered.bam$",
        full.names = TRUE
    )

    if (length(aligned_files) < 2) {
        stop(paste(
            "Non aligned BAM files with the extension _filtered.bam in",
            align_dir
        ))
    }

    alignedR1_files <- aligned_files[grep("_R1", aligned_files)]
    alignedR2_files <- aligned_files[grep("_R2", aligned_files)]

    # Load digested genome
    # Pattern matching breaks for ensembl based genomes (1.rda -> 11.rda)
    # Assume that viewpoint can only correspond to 1 chromosome anyway, and fail if the 
    file <- file.path(digested_genome, paste0(as.character(GenomicRanges::seqnames(pos_viewpoint)), ".rda"))
    if (!file.exists(file)) stop(paste("No digested genome file found at: ", file))

    #file <- list.files(digested_genome,
    #    pattern = paste0(as.character(GenomicRanges::seqnames(pos_viewpoint)), ".rda"),
    #    full.names = TRUE
    #)
    #if (length(file) == 0) stop(paste("No digested genome file"))

    digested_genome_gr <- NULL # Avoid no visible binding for global variable
    load(file)


    nll <- lapply(
        seq_len(length(alignedR1_files)),
        function(i) {
            .singleCounterUMI4C(
                filtered_bam_R1 = alignedR1_files[i],
                filtered_bam_R2 = alignedR2_files[i],
                digested_genome_gr = digested_genome_gr,
                pos_viewpoint = pos_viewpoint,
                res_enz = res_enz,
                count_dir = count_dir,
                filter_bp = filter_bp
            )
        }
    )
}


#' Count UMIs for a given bam file.
#' @param filtered_bam_R1 R1 bam file.
#' @param filtered_bam_R2 R2 bam file.
#' @param digested_genome_gr GRanges object containing the coordinates for the
#' digested genome.
#' @param pos_viewpoint Vector consist of chromosome, start and end position of
#' the viewpoint.
#' @param count_dir Counter directory.
#' @return Creates a tab-delimited file in \code{wk_dir/count} named
#' "\code{basename(fastq) _counts.tsv}", containing the
#' coordinates for the viewpoint fragment, contact fragment and the number of
#' UMIs detected in the ligation.
#' @inheritParams contactsUMI4C
.singleCounterUMI4C <- function(filtered_bam_R1,
    filtered_bam_R2,
    digested_genome_gr,
    pos_viewpoint,
    res_enz,
    count_dir,
    filter_bp = 10e6) {
    viewpoint_filter <- GenomicRanges::resize(pos_viewpoint,
        width = filter_bp * 2,
        fix = "center"
    )

    # load digest genome, filter and transform to granges
    digested_genome_gr <- IRanges::subsetByOverlaps(
        digested_genome_gr,
        viewpoint_filter
    )

    # read bam and transform to a granges
    bam_R1_gr <- GenomicAlignments::readGAlignments(filtered_bam_R1,
        use.names = TRUE,
        param = Rsamtools::ScanBamParam(which = viewpoint_filter)
    )
    bam_R2_gr <- GenomicAlignments::readGAlignments(filtered_bam_R2,
        use.names = TRUE,
        param = Rsamtools::ScanBamParam(which = viewpoint_filter)
    )

    mcols(bam_R1_gr)$header <- names(bam_R1_gr)
    mcols(bam_R1_gr)$umi <- vapply(
        strsplit(names(bam_R1_gr), ":"),
        function(x) x[1],
        FUN.VALUE = character(1)
    )
    mcols(bam_R1_gr)$readID <- vapply(
        strsplit(names(bam_R1_gr), ":"),
        function(x) paste0(x[2], "_", x[3]),
        FUN.VALUE = character(1)
    )

    mcols(bam_R2_gr)$header <- names(bam_R2_gr)
    mcols(bam_R2_gr)$umi <- vapply(
        strsplit(names(bam_R2_gr), ":"),
        function(x) x[1],
        FUN.VALUE = character(1)
    )
    mcols(bam_R2_gr)$readID <- vapply(
        strsplit(names(bam_R2_gr), ":"),
        function(x) paste0(x[2], "_", x[3]),
        FUN.VALUE = character(1)
    )


    gr1 <- GenomicAlignments::granges(bam_R1_gr, use.mcols = TRUE, use.names = FALSE)
    gr2 <- GenomicAlignments::granges(bam_R2_gr, use.mcols = TRUE, use.names = FALSE)

    gr <- c(gr1, gr2)

    gr_sp <- GenomicRanges::GRangesList(GenomicRanges::split(gr, gr$readID))

    n <- elementNROWS(gr_sp)

    ## Regions containing only the viewpoint
    uni <- unlist(gr_sp[n == 1])
    uni <- regioneR::joinRegions(uni)

    ligations <- gr_sp[n > 1]
    ligations_noview <- IRanges::subsetByOverlaps(unlist(ligations), uni,
        invert = TRUE
    )
    lig_noview_unl <- GenomicRanges::GRangesList(GenomicRanges::split(
        ligations_noview,
        ligations_noview$readID
    ))

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
        minoverlap = nchar(res_enz) + 1
    )

    final <- final_ligations[queryHits(hits)]
    final$fragID <- as.character(mcols(digested_genome_gr)[subjectHits(hits), 1])
    table <- as.data.frame(table(final$fragID))

    umis_df <- data.frame(digested_genome_gr)[, c(1, 2, 3, 6)]
    umis_df <- suppressWarnings(dplyr::left_join(umis_df, table, by = c(id = "Var1")))

    viewpoint <- data.frame(subsetByOverlaps(digested_genome_gr, pos_viewpoint))[, c(1, 2, 3)]

    final_umis <- cbind(
        viewpoint[rep(1, nrow(umis_df)), ],
        umis_df[, -4]
    )
    colnames(final_umis) <- c(
        "chr_bait", "start_bait", 'end_bait',
        "chr_contact", "start_contact", 'end_contact',
        "UMIs"
    )
    final_umis$UMIs[is.na(final_umis$UMIs)] <- 0

    file_name <- strsplit(basename(filtered_bam_R1), "_R1")[[1]][1]
    counts_file <- file.path(count_dir, paste0(file_name, "_counts.tsv"))

    # gz file directly to avoid errors due to disk latency
    gzf <- gzfile(paste0(counts_file, '.gz'), "w")
    utils::write.table(
    x = final_umis,
    file = gzf,
    row.names = FALSE,
    quote = FALSE,
    sep = "\t"
    )
    close(gzf)


    # R.utils::gzip(counts_file, overwrite = TRUE)

    message(
        paste0("[", Sys.time(), "] "),
        "Finished sample ", file_name
    )
}
