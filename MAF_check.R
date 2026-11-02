args <- commandArgs(trailingOnly = TRUE)
location <- args[1]

maf_freq <- read.table(file=paste0(location, "MAF_check.frq"), header =TRUE, as.is=T)
pdf("MAF_distribution.pdf")
hist(maf_freq[,5],main = "MAF distribution", xlab = "MAF")
dev.off()