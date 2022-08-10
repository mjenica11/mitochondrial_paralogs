# Convert salmon quant.sfs to count matrix; use these counts to generate lms

library(tximport)

# Read in salmon quants
directories <- Sys.glob("SRR74267*")

files <- file.path("/scratch/mjpete11/sweet_data", directories, "quant.sf")

names(files) <- directories

# Some of the quant files are missing so drop these files for the moment
files <- files[-c(4,5,11,14)]

# Create tx2gene object
library(GenomicFeatures)
txdb <- makeTxDbFromGFF("/scratch/mjpete11/sweet_data/gencode.v41.annotation.gtf", format=c("auto", "gff3", "gtf"))
k <- keys(txdb, keytype = "TXNAME")
tx2gene <- select(txdb, k, "GENEID", "TXNAME")

txi.salmon <- tximport(files, type="salmon", tx2gene=tx2gene)
