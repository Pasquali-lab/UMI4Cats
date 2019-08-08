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

colData <- data.frame(sampleID = gsub(".txt", "", files),
                      replicate = unlist(lapply(strsplit(files, "_"), function(x) x[3])),
                      condition = gsub(".txt", "", unlist(lapply(strsplit(files, "_"), function(x) x[4]))),
                      file = paste0("data-raw/", files),
                      stringsAsFactors=F)

UMI4C <- UMI4C(colData,
            viewpoint_name="yes")
