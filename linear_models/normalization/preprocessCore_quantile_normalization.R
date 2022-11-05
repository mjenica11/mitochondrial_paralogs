# Try quantile normalization with a different package

# Libraries
library(data.table)
library(tidyverse)
library(reshape)
library(preprocessCore)

#_______________________________________________________________________________ 
# Read in combat_seq batch adjusted matrix 
#_______________________________________________________________________________ 
# Read in GTEx counts
counts <- fread("/scratch/mjpete11/linear_models/data/combat_seq_filtered.csv", sep=",")

# Drop the index column
counts$V1 <- NULL

#_______________________________________________________________________________ 
# Perform quantile normalization 
#_______________________________________________________________________________ 
#perform quantile normalization
mat  <- as.matrix(counts[,2:ncol(counts)])
df_norm <- as.data.frame(normalize.quantiles(mat))

df_norm[1:5,1:5]
counts[1:5,2:6]

dim(df_norm)==dim(counts[,2:ncol(counts)]) # TRUE TRUE

# Add sample names column
colnames(df_norm) <- colnames(counts[,2:ncol(counts)])

# Add gene name column
gene_names <- counts$"Description" 

# Add gene name column 
df_norm$"Description" <- gene_names

# Move the gene name column to the front
df_norm2 <- df_norm %>% select("Description", everything())

df_norm2[1:5,1:5]
counts[1:5,1:5]

# Write to file
write.csv(df_norm2, "/scratch/mjpete11/linear_models/data/preprocessCore_quantile_normalized_counts.csv")
