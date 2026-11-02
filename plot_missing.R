args <- commandArgs(trailingOnly = TRUE)
location <- args[1]
filename <- args[2]

# read data into R 
indmiss<-read.table(file=paste0(location, filename, ".imiss"), header=TRUE)
snpmiss<-read.table(file=paste0(location, filename, ".lmiss"), header=TRUE)


pdf("histimiss.pdf") #indicates pdf format and gives title to file
hist(indmiss[,6],main="Histogram individual missingness") #selects column 6, names header of file

pdf("histlmiss.pdf") 
hist(snpmiss[,5],main="Histogram SNP missingness")  

dev.off() # shuts down the current device


