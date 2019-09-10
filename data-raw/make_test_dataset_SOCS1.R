library(umi4cPackage)
p4cLoadConfFiles("~/Projects/CYT_hg19/data/UMI4C/sep/conf/")

t <- gtrack.ls()
sel <- t[grepl("SOCS", t) & !grepl("_m_", t)]

for (i in sel) {
  t <- p4cNewProfile(i,
                     scope_5=5e6, scope_3=5e6)
  df <- as.data.frame(t$dgram[,1:2])
  colnames(df) <- c("pos_contact", "UMIs")
  df$chr_bait <- paste0("chr", t$bait$chrom)
  df$chr_contact <- df$chr_bait
  df$pos_bait <- as.numeric(t$bait$start)
  df <- df[,c("chr_bait", "pos_bait", "chr_contact", "pos_contact", "UMIs")]

  write.table(df, file=paste0("data-raw/", i, ".txt"), row.names=F, sep="\t", quote=F)
}


#### Testing object class UMI4C

files <- list.files("data-raw",
                    pattern="umi4C_SOCS1*",
                    full.names=F)
# files <- files[-c(2:3)]

colData <- data.frame(sampleID = gsub(".txt", "", files),
                      replicate = unlist(lapply(strsplit(files, "_"), function(x) x[3])),
                      condition = gsub(".txt", "", unlist(lapply(strsplit(files, "_"), function(x) x[4]))),
                      file = paste0("data-raw/", files),
                      stringsAsFactors=F)

UMI4C <- UMI4C(colData,
            viewpoint_name="SOCS1")

UMI4C <- processUMI4C(UMI4C,
                      min_win_cov=NULL)

UMI4C <- normalizeUMI4C(UMI4C,
                        ref_umi4c="smaller")

ggplot(metadata(UMI4C)$trend_norm) +
  geom_line(aes(coord, trend,
                group=interaction(group, sample),
                color=sample))


dgram <- as.data.frame(metadata(UMI4C)$norm_dgram$umi4C_SOCS1_HI22_ctrl)
dgram$start <- start(rowRanges(UMI4C))
dgram$end <- (dgram$start[c(2:nrow(dgram), nrow(dgram))] -
                dgram$start) + dgram$start

dgram.l <- reshape2::melt(dgram,
                          id.vars=(ncol(dgram)-1):ncol(dgram),
                          value.vars=1:(ncol(dgram)-2),
                          variable.name="scales",
                          value.name="value")
dgram.l$scales <- as.numeric(dgram.l$scales)

## Log2 + normalize by maximum value in view
norm_fact <- max(dgram.l$value, na.rm=T)
dgram.l$value_norm <- log2(dgram.l$value+1) - log2(norm_fact+1)

## Plot dgram
# domainogram <-
  ggplot2::ggplot(dgram.l) +
  ggplot2::geom_rect(ggplot2::aes(xmin=start, xmax=end,
                                  ymin=rev(scales), ymax=rev(scales)+1,
                                  fill=value_norm)) +
  viridis::scale_fill_viridis(option = "B",
                              na.value=NA,
                              labels=function(x) paste0(round(100*2^x), "%"),
                              name="Contacts/Maximum",
                              breaks=scales::pretty_breaks(n=4),
                              guide = ggplot2::guide_colorbar(direction = "horizontal",
                                                              title.position="top",
                                                              barwidth=8))
