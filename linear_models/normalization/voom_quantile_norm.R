# Try quantile normalization and logCPM transformation with limma voom

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

# Generate the design matrix
# Samples are rows and columns are covariates of interest (heart and liver)
# Use this file to create the design matrix
organs <- read.csv("/scratch/mjpete11/linear_models/data/organs_combat_seq_limma.csv", sep=",")
organs$value <- NULL
organs$SUBJID <- NULL
organs$gene <- NULL
organs$SMRIN <- NULL
organs$SMTSISCH <- NULL
organs2 <- organs[!duplicated(organs),]

# Drop samples that are not from paired heart or liver
keeps <- as.vector(organs2$SAMPID)
counts2 <- subset(counts, select=keeps) 

ncol(counts2)==nrow(organs2)

# Reshape df so it can be in the correct form to make a design matrix
organs3 <- table(organs2$SAMPID, organs2$organ)

# Reorder the samples (columns) to match the order of the samples
# in the design matrix
sorted_organs2 <- organs2[order(match(organs2$SAMPID,row.names(organs3))),]

# Double check that the order of sample IDs in sorted_organs2 is equal to
# the order of the row names in organs3
identical(row.names(organs3), sorted_organs2$SAMPID) # TRUE

# Convert table to matrix for use with voom function
rownames(organs3) <- c() # drop row names
mat <- as.matrix(organs3)

#_______________________________________________________________________________ 
# Perform quantile normalization 
#_______________________________________________________________________________ 
# Perform quantile normalization with voom
# Transpose since the function applies qnrom by column
# but the untransposed matrix is genes (rows) by samples (columns)
mat_counts <- as.matrix(counts2)
ncol(mat_counts)==nrow(organs3) # TRUE

# Apply voom normalization
#df_norm <- as.data.frame(voom(counts=mat_counts[1:20,1:20], design=mat[1:20,], normalize.method="quantile", save.plot=TRUE)) # TEST
df_norm <- voom(counts=mat_counts, design=mat, normalize.method="quantile", save.plot=TRUE)

df_norm[1:5,1:5]
counts2[1:5,2:6]

#dim(df_norm)==dim(counts[,2:ncol(counts)]) # TRUE TRUE
dim(df_norm[["weights"]])==dim(counts2) # TRUE TRUE

# Save the logCPM and voom normalized counts into a separate df
dat <- as.data.frame(df_norm[["weights"]])

# Add column names back
colnames(dat) <- colnames(df_norm$E) 

# Add gene name column
gene_names <- counts$"Description"

# Add gene name column 
dat$"Description" <- gene_names

# Move the gene name column to the front
dat2 <- dat %>% select("Description", everything())

dat2[1:5,1:5]
counts2[1:5,1:5]

# Double check that the column names of the original counts and voom
# transformed counts are in the same order
identical(colnames(counts2), colnames(dat2[,2:ncol(dat2)])) # TRUE

# Write to file
write.csv(dat2, "/scratch/mjpete11/linear_models/data/voom_quantile_normalized_counts.csv")
