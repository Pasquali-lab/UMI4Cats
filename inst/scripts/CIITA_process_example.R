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
  out_path = "digested_genome/"
)

## 2) Process UMI-4C fastq files -----------------------------------------------
raw_dir <- file.path(path, "CIITA", "fastq")

contactsUMI4C(
  fastq_dir = raw_dir,
  wk_dir = "CIITA",
  bait_seq = "GGACAAGCTCCCTGCAACTCA",
  bait_pad = "GGACTTGCA",
  res_enz = "GATC",
  cut_pos = 0,
  digested_genome = hg19_dpnii,
  bowtie_index = file.path(path, "ref_genome", "ucsc.hg19.chr16"),
  ref_gen = BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19,
  threads = 10
)
