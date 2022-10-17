# Try quantile normalization with a different package

# Libraries
library(data.table)
library(tidyverse)
library(reshape)
library(sva)
library(preprocessCore)

#_______________________________________________________________________________ 
# Make matrix of genes, samples, counts, organs, and batch effects
#_______________________________________________________________________________ 
# Read in GTEx manifest
manifest <- read.csv("/scratch/mjpete11/linear_models/data/sample.tsv", header=TRUE, sep = "\t")

# Make dataframe with sample id, tissue type
# All of the rin number and ischemic time values were missing...
df1 <- data.frame(manifest$"dbgap_sample_id", manifest$"tissue_type")

# Remove all rows with NA in either columns
df2 <- df1[complete.cases(df1), ]

# Read in GTEx counts
counts <- fread("/scratch/mjpete11/linear_models/data/combat_batch_adjusted_counts.csv", sep=",")

# Drop the index column
counts$V1 <- NULL

#_______________________________________________________________________________ 
# Perform quantile normalization 
#_______________________________________________________________________________ 
#perform quantile normalization
df_norm <- as.data.frame(normalize.quantiles(as.matrix(counts)))

df_norm[1:5,1:5]
counts[1:5,1:5]

dim(df_norm)==dim(counts) # TRUE TRUE

# Add sample names column
colnames(df_norm) <- colnames(counts)

# Add gene name column
gene_names <- fread("/scratch/mjpete11/linear_models/data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct", sep="\t",
					   	select=c("Description"))
df_norm$"Description" <- gene_names 

# Move the gene name column to the front
df_norm2 <- df_norm %>% select("Description", everything())

df_norm2[1:5,1:5]
counts[1:5,1:5]

# Write to file
write.csv(df_norm2, "/scratch/mjpete11/linear_models/data/preprocessCore_quantile_normalized_counts.csv")
