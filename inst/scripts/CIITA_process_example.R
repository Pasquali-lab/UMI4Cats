##################################
# PROCESS CIITA EXAMPLE UMI4Cats #
##################################
# This script describes how files contained in the following folders were
# generated.
#    - inst/extdata/CIITA/count
#    - inst/extdata/CIITA/logs
#
# Briefly, CIITA processed files were generated using the UMI4Cats package
# and the sample dataset that can be downloaded using downloadUMI4CexampleData().
# For more information on the origin of fastq files, please check:
#   - inst/scripts/CIITA_generate_reduced_fastq.R

## 0) Load package & download example data -------------------------------------
library(UMI4Cats)

path <- downloadUMI4CexampleData()

## 1) Generate Digested genome -------------------------------------------------
hg19_dpnii <- digestGenome(
  cut_pos = 0,
  res_enz = "GATC",
  name_RE = "DpnII",
  ref_gen = BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19,
  sel_chr = "chr16",
  out_path = "digested_genome/"
)

## 2) Process UMI-4C fastq files -----------------------------------------------
raw_dir <- file.path(path, "CIITA", "fastq")
wk_dir <- here::here("inst", "extdata", "CIITA")

contactsUMI4C(
  fastq_dir = raw_dir,
  wk_dir = wk_dir,
  bait_seq = "GGACAAGCTCCCTGCAACTCA",
  bait_pad = "GGACTTGCA",
  res_enz = "GATC",
  cut_pos = 0,
  digested_genome = hg19_dpnii,
  bowtie_index = file.path(path, "ref_genome", "ucsc.hg19.chr16"),
  ref_gen = BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19,
  filter_bp = 1e6,
  threads = 10
)

## 3) Create UMI4C object -----------------------------------------------------
# Load sample processed file paths
files <- list.files(file.path(wk_dir, "count"),
                    pattern = "*_counts.tsv.gz",
                    full.names = TRUE
)

# Create colData including all relevant information
colData <- data.frame(
  sampleID = gsub("_counts.tsv.gz", "", basename(files)),
  file = files,
  stringsAsFactors = FALSE
)

library(tidyr)
colData <- colData %>%
  separate(sampleID,
           into = c("condition", "replicate", "viewpoint"),
           remove = FALSE
  )

# Load UMI-4C data and generate sample UMI4C object
ex_ciita_umi4c <- makeUMI4C(
  colData = colData,
  viewpoint_name = "CIITA",
  bait_expansion = 3e5,
  grouping = NULL
)

usethis::use_data(ex_ciita_umi4c, overwrite=TRUE)
