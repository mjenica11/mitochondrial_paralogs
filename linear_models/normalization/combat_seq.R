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

# Expression filter: Keep a gene if it has >5 counts in at least 1 samples 
expression_filter <- function(DF, thresh){
				DAT <- DF[rowSums(DF>thresh)>=1, ]
				return(DAT)
}

# Apply function
counts3 <- expression_filter(DF=counts2, thresh=5)

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
setdiff(SLC, counts3$'Description') 

# How many genes are left?
nrow(counts3) # 

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
