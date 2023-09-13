# Perform combat normalization on count data

# Libraries
library(data.table)
library(tidyverse)
library(reshape)
library(sva)
library(dplyr)

# Read in GTEx counts
counts <- fread("/scratch/mjpete11/linear_models/data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct", sep="\t")

# Read in metadata with batch variables
manifest <- read.csv("/scratch/mjpete11/linear_models/data/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt", header=TRUE, sep = "\t")

# Drop samples from manifest that are missing in the count data
rownames(manifest) <- manifest$SAMPID
manifest <- manifest[(manifest$SAMPID %in% colnames(counts)[3:ncol(counts)]),]
samples <- manifest$SAMPID

# Only keep samples that are present in the manifest 
counts2 <- counts[,..samples]

##### TEST ####
#counts2 <- counts2[1:10,1:10]
##### TEST ####

# Replace 0s with 1s 
counts2[counts2 == 0] <- 1

# Log transform the counts but skip the first column
counts2 <- as.data.frame(log10(as.matrix(counts2)))

# Check if any negative values
any(counts2)<0 # FALSE

##### TEST ####
#counts2$"Description" <- c("Hello","Hello","Hello","Hello","Hello","Hello","Hello","Hello","Hello","Hello")
##### TEST ####

# Add gene names column back on
counts2$"Description" <- counts$"Description"

# Move the gene names column to the front
counts2 <- counts2 %>% select("Description", everything())

# Median filter: Function to filter rows that do not have a median >= condition
median_filter <- function(DF, thresh){
				DAT <- DF[which(apply(DF[,-c(1)],1,median) > thresh), ]
				return(DAT)
}

# Filter rows with a median of < 0 log(counts)
counts3 <- median_filter(DF=counts2, thresh=0)

# Drop samples from manifest that are missing in the count data
mani <- manifest[(manifest$SAMPID %in% colnames(counts3)[2:ncol(counts3)]),]

# Make a list of known batch variables; genotype or expression batch ID 
batch <- mani$SMGEBTCH

# Create a model matrix for the adjustment variables
# Only added and intercept term, since I am only adjusting for batch
modcombat <- model.matrix(~1, data=mani)

# Apply ComBat to the data, using parametric empirical Bayesian adjustment
combat_edata <- ComBat(dat=counts3[,-c(1)], batch=batch, mod=modcombat, par.prior=TRUE)

# Append gene names back 
combat_edata <- as.data.frame(combat_edata)
combat_edata$"Description" <- counts3$"Description"

# Move the gene names column to the front
combat_edata <- combat_edata %>% select("Description", everything())

# Write to file
write.csv(combat_edata, "/scratch/mjpete11/linear_models/data/combat_batch_adjusted_counts_zero_filtered.csv")
