% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/digestGenome.R
\name{digestGenome}
\alias{digestGenome}
\title{Digest reference genome}
\usage{
digestGenome(cut_seq_5p, cut_seq_3p, name_RE, refgen,
  out_path = "digested_genome/")
}
\arguments{
\item{cut_seq_5p}{5' restriction enzyme cut sequence. For example, for DpnII |GATC \code{cut_seq_5p=""}.}

\item{cut_seq_3p}{3' restriction enzyme cut sequence. For example, for DpnII |GATC \code{cut_seq_3p="GATC"}.}

\item{name_RE}{Restriction enzyme name .}

\item{refgen}{A BSgenome object of the reference genome.}

\item{out_path}{Output path where to save the genomic track. The default is a directory named \code{digested_genome/}
created in your working directory.}
}
\value{
Creates a tab-delimited file in \code{out_path} named "\code{refgen_name_RE.tsv}", containing the
coordinates for all fragments generated by digesting the genome using the given cutting sequences.
}
\description{
Performs an \emph{in silico} digestion of a given reference genome using a given restriction enzyme sequence.
}
\examples{
\dontrun{
library(BSgenome.Hsapiens.UCSC.hg19)
refgen <- BSgenome.Hsapiens.UCSC.hg19

hg19_dpnii <- digestGenome(cut_seq_5p="",
                           cut_seq_3p="GATC",
                           name_RE="DpnII",
                           refgen=refgen,
                           out_path="digested_genome/")
}
}