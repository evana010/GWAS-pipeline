args <- commandArgs(trailingOnly = TRUE)
location <- args[1]

hwe<-read.table(paste0(location, "hardy.hwe"), header=TRUE)
pdf("histhwe.pdf")
hist(hwe[,9],main="Histogram HWE")
dev.off()

