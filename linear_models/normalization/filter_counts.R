# Perform combat normalization on count data

# Libraries
library(data.table)
library(tidyverse)
library(reshape)
library(sva)
library(dplyr)

# Read in GTEx counts
counts <- fread("/scratch/mjpete11/linear_models/data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct", sep="\t")

dim(counts) # 56200 17384

# Read in metadata with batch variables
manifest <- read.csv("/scratch/mjpete11/linear_models/data/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt", header=TRUE, sep = "\t")

# Drop samples from manifest that are missing in the count data
rownames(manifest) <- manifest$SAMPID
manifest <- manifest[(manifest$SAMPID %in% colnames(counts)[3:ncol(counts)]),]

# ComBat-seq does not support 1 sample per batch
# Remove samples (rows) that were sequenced alone from the manifest
mani <- manifest[manifest$SMGEBTCH %in% manifest$SMGEBTCH[duplicated(manifest$SMGEBTCH)],]

# Only keep samples that are present in the manifest which were sequenced 
# with at least one other sample
samples <- mani$SAMPID
counts2 <- counts[,..samples]

# Append gene names back 
counts2$"Description" <- counts$"Description"

# Move the gene names column to the front
counts2 <- counts2 %>% select("Description", everything())
dim(counts2) # 56200 17312

# Expression filter: Keep a gene if it has >x counts in at least 1 sample 
# Skip the first column because that is just the gene name column
expression_filter <- function(DF, thresh){
				DAT <- DF[rowSums(DF[,2:ncol(DF)]>=thresh, na.rm=TRUE)>0, ]
				return(DAT)
}

# Apply function
counts3 <- expression_filter(DF=counts2, thresh=5)
dim(counts3) # 51259 17312 

# Use another function to test if it is working...
tmp <- counts2[rowSums(counts2[,2:ncol(counts2)]>=5, na.rm=TRUE)>0,] 
dim(tmp) # 51259 17312

################################ test code ####################################
# Create a sample dataframe
dat <- data.frame(a = c("X", "Y", "Z", "W"),
				 b = c(1, 2, 3, 4),
				 c = c(5, 4, 3, 2),
				 d = c(9, NA, 11, 12))

# Drop rows with no value greater than or equal to 5 in columns b, c, and d
dat <- dat[rowSums(dat[, 2:ncol(dat)] >= 5, na.rm = TRUE) > 0, ]
dim(dat) #  3 4
# Drop rows with no value greater than or equal to 5
dat1 <- dat[rowSums(dat[, 2:ncol(dat)] >= 5, na.rm = TRUE) > 0, ]
dim(dat1) #  3 4
dat2 <- expression_filter(DF=dat, thresh=5)
dim(dat2) # 3 4
#...Both functions are definitely working...
# Diagnosed the problem! Needed to explicitely skip the gene name column
################################ test code ####################################

# How many genes are left if you set the threshold to 10?
counts_10 <- expression_filter(DF=counts2, thresh=10)
dim(counts_10) # 47465 17312
test0 <- counts2[rowSums(counts2[,2:ncol(counts2)]>=10, na.rm=TRUE)>0,] # Still got the same result
dim(test0) # 47465 17313.. results can be replicated

# Write filtered counts to file
write.csv(counts3, "/scratch/mjpete11/linear_models/data/filtered_counts.csv")

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
setdiff(SLC, counts3$'Description')  # No
