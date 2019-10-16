devtools::load_all()

files <- list.files("data-raw",
                    pattern="umi4C_SOCS1*",
                    full.names=F)

colData <- data.frame(sampleID = gsub(".txt", "", files),
                      replicate = unlist(lapply(strsplit(files, "_"), function(x) x[3])),
                      condition = gsub(".txt", "", unlist(lapply(strsplit(files, "_"), function(x) x[4]))),
                      file = paste0("data-raw/", files),
                      stringsAsFactors=F)

umi <- makeUMI4C(colData=colData,
                 viewpoint_name="SOCS1",
                 min_win_factor=0.02,
                 normalized=FALSE)

umi_norm <- makeUMI4C(colData=colData,
                      viewpoint_name="SOCS1",
                      min_win_factor=0.02,
                      normalized=TRUE)

plotUMI4C(umi,
          grouping="condition")
plotUMI4C(umi_norm,
          trend_grouping="condition",
          ylim=c(0,5),
          xlim=c(11.0e6, 11.55e6))

p2 <- plotUMI4C(umi_norm,
                grouping="condition",
                ylim=c(0,5))

cowplot::plot_grid(p2, p1, ncol=1, align="v")
