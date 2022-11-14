# Try quantile normalization with a different package

# Libraries
library(data.table)
library(tidyverse)
library(reshape)
library(limma)

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
# Perform quantile normalization
# Transpose since the function applies qnrom by column
# but the matrix is genes (rows) by samples (columns)
mat  <- as.matrix(t(counts[,2:ncol(counts)]))
df_norm <- as.data.frame(normalizeQuantiles(mat))

df_norm[1:5,1:5]
counts[1:5,2:6]

dim(df_norm)==dim(t(counts[,2:ncol(counts)])) # TRUE TRUE

# Add sample names column
df_norm2 <- as.data.frame(t(df_norm))
dim(df_norm2)==dim(counts[,2:ncol(counts)]) # TRUE TRUE
colnames(df_norm2) <- colnames(counts[,2:ncol(counts)])
df_norm2[1:5,1:5]

# Add gene name column
gene_names <- counts$"Description"

# Add gene name column 
df_norm2$"Description" <- gene_names

# Move the gene name column to the front
df_norm3 <- df_norm2 %>% select("Description", everything())

df_norm3[1:5,1:5]
counts[1:5,1:5]

# Write to file
write.csv(df_norm3, "/scratch/mjpete11/linear_models/data/limma_quantile_normalized_counts.csv")
