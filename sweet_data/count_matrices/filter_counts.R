# Filter counts after combining into one matrix

# Libraries
library(data.table)
library(tidyverse)
library(reshape)
library(sva)
library(dplyr)

# Set the path
path <- paste0(getwd(),"/")
path

# Read in count matrix
counts <- read.csv(paste0(path, "/combined_filtered_matrix.csv")) 
counts[1:5,1:5]
counts$X <- NULL
dim(counts) # 23346 65 

# Read in metadata with batch variables
manifest <- read.csv(paste0(path, "/metadata.csv"), header=TRUE, sep = ",")
dim(manifest) # 64 30
head(manifest)

# Samples to drop from metadata bc Salmon couldn't process it: 
# RUN: SRR7426822
# BioSample: SAMN09484097 
manifest <- manifest[!manifest$BioSample=="SAMN09484097",]
dim(manifest) # 63 30

# Move the gene names column to the front
counts <- counts %>% select("Hugo_ID", everything())
counts[1:5,1:5]

# Expression filter: Keep a gene if it has >x counts in at least 10 samples
# Skip the first column because that is just the gene name column
expression_filter <- function(DF, thresh){
				DAT <- DF[rowSums(DF[,2:ncol(DF)]>=thresh, na.rm=TRUE)>10, ]
				return(DAT)
}

# Apply function
counts_2 <- expression_filter(DF=counts, thresh=5)
dim(counts_2) # 22540 64 
dim(counts) # 23346 64 ; dropped 806 genes 
counts_2[1:5,1:5]
counts[1:5,1:5]

# Use another function to test if it is working...
tmp <- counts[rowSums(counts[,2:ncol(counts)]>=5, na.rm=TRUE)>0,] 
dim(tmp) # 23346 14 

# Write filtered counts to file
write.csv(counts_2, paste0(path, "/total_filtered_counts.csv"))

# Did any SLCs drop out?
SLC <- c("SLC25A1", "SLC25A2", "SLC25A3", "SLC25A4", "SLC25A5", "SLC25A6",
		 "UCP1", "UCP2", "UCP3", "SLC25A10", "SLC25A11", "SLC25A12", 
		 "SLC25A13", "SLC25A14", "SLC25A15", "SLC25A16", "SLC25A17", 
		 "SLC25A18", "SLC25A19", "SLC25A20", "SLC25A21", "SLC25A22", 
		 "SLC25A23", "SLC25A24", "SLC25A25", "SLC25A26", "SLC25A27", 
		 "SLC25A28", "SLC25A29", "SLC25A30", "SLC25A31", "SLC25A32", 
		 "SLC25A33", "SLC25A34", "SLC25A35", "SLC25A36", "SLC25A37", 
		 "SLC25A38", "SLC25A39", "SLC25A40", "SLC25A41", "SLC25A42", 
		 "SLC25A43", "SLC25A44", "SLC25A45", "SLC25A46", "SLC25A47",
		 "SLC25A48", "MTCH1", "MTCH2", "SLC25A51", "SLC25A52", "SLC25A53")

# Are any genes missing?
setdiff(SLC, counts_2$'Hugo_ID')  # "SLC25A2" "UCP1" "SLC25A31" "SLC25A47" 
