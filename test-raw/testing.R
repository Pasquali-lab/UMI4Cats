# library(UMI4Cats)
devtools::load_all()

files <- list.files("data-raw",
                    pattern="umi4C_SOCS1*",
                    full.names=F)

# files <- files[c(1,3)]
# files <- files[1]

colData <- data.frame(sampleID = gsub(".txt", "", files),
                      replicate = unlist(lapply(strsplit(files, "_"), function(x) x[3])),
                      condition = gsub(".txt", "", unlist(lapply(strsplit(files, "_"), function(x) x[4]))),
                      file = paste0("data-raw/", files),
                      stringsAsFactors=F)


umi4c <- makeUMI4C(colData=colData,
                      viewpoint_name="SOCS1",
                      min_win_factor=0.02,
                      normalized=TRUE)

plotUMI4C(umi4c,
          grouping="condition", #c("condition", "replicate"),
          dgram_function="quotient",
          dgram_plot=F,
          protein_coding=F,
          ylim=c(0,5),
          xlim=c(11.0e6, 11.45e6),
          font_size=12)


query_regions <- GenomicRanges::GRanges(c("chr16:11247500-11252500",
                                          "chr16:11397500-11402500",
                                          "chr16:11300000-11305000"))

umi4c <- fisherUMI4C(umi4c,
                          query_regions=query_regions)

umi4c <- fisherUMI4C(umi4c,
                     padj_method = "BH")
plotUMI4C(umi4c, dgram_plot=T)

plotUMI4C(umi4c,
          grouping="condition", #c("condition", "replicate"),
          dgram_function="quotient",
          dgram_plot=F,
          protein_coding=F,
          ylim=c(0,5),
          xlim=c(11.0e6, 11.45e6),
          font_size=12)

umi4c <- deseq2UMI4C(umi4c)

plotUMI4C(umi4c,
          grouping="condition", # c("condition", "replicate"),
          dgram_function="quotient",
          dgram_plot=T,
          protein_coding=F,
          ylim=c(0,5),
          xlim=c(11.0e6, 11.45e6),
          font_size=12)
