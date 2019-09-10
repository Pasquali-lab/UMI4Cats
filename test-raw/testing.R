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
                 min_win_factor=0.02)

metadata(umi)
dgram(umi)

plot(assays(umi)$geo_coords[,1], assays(umi)$trend[,1], type="l")
## Make norm mat
# dgram * norm_mat
# recalculate trend
scaled <- calculateAdaptativeTrend(umi, scaled=T)
