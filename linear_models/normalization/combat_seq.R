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

# Median filter: Drop a gene if no individual has <5 counts
median_filter <- function(DF, thresh){
				DAT <- DF[!rowSums(DF<5), ]
				return(DAT)
}

# Filter rows with a median of < 0 log(counts)
# Skip first column when applying function since it contains the gene names
counts3 <- median_filter(DF=counts2, thresh=5)

# How many genes are left?
nrow(counts3)

# None of the batch IDs should be unique
length(unique(mani$SMGEBTCH))

# Make a list of known batch variables; genotype or expression batch ID 
batch <- mani$SMGEBTCH

# Create a model matrix for the adjustment variables
# Only added and intercept term, since I am only adjusting for batch
modcombat <- model.matrix(~1, data=mani)

# Apply ComBat_seq to the data, using parametric empirical Bayesian adjustment
combat_edata <- ComBat_seq(counts=as.matrix(counts3[,-c(1)]), batch=batch)

# Append gene names back 
combat <- as.data.frame(combat_edata)
combat$"Description" <- counts3$"Description"

# Move the gene names column to the front
combat <- combat %>% select("Description", everything())

# Write to file
write.csv(combat, "/scratch/mjpete11/linear_models/data/combat_seq_filtered.csv")
