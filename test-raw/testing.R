library(UMI4Cats)
# devtools::load_all()

files <- list.files("data-raw",
                    pattern="umi4C_SOCS1*",
                    full.names=F)

# files <- files[1]

colData <- data.frame(sampleID = gsub(".txt", "", files),
                      replicate = unlist(lapply(strsplit(files, "_"), function(x) x[3])),
                      condition = gsub(".txt", "", unlist(lapply(strsplit(files, "_"), function(x) x[4]))),
                      file = paste0("data-raw/", files),
                      stringsAsFactors=F)


umi_norm <- makeUMI4C(colData=colData,
                      viewpoint_name="SOCS1",
                      min_win_factor=0.02,
                      normalized=TRUE)

plotUMI4C(umi_norm,
          grouping="condition", #c("condition", "replicate"),
          dgram_function="quotient",
          dgram_plot=F,
          protein_coding=F,
          ylim=c(0,5),
          xlim=c(11.0e6, 11.45e6),
          font_size=12)


