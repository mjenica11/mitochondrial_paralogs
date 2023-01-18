# Try quantile normalization and logCPM transformation with limma voom

# Libraries
library(data.table)
library(tidyverse)
library(reshape)
library(limma)
library(plyr)
library(stringr)
library(ggplot2)
library(ggpubr)

# Read in the organs metadataframe created in the batch_voom_qnorm.R script
# Use this file to create the design matrix
organs <- read.csv("/scratch/mjpete11/linear_models/data/filtered_counts_organ_metadata.csv", sep=",")

# Read in counts
counts <- fread("/scratch/mjpete11/linear_models/data/filtered_counts.csv", sep=",")

# Drop the index column
counts$V1 <- NULL

# Drop samples from the counts df that are not from paired heart or liver
keeps <- as.vector(organs$SAMPID)
counts2 <- subset(counts, select=keeps) 

# Are there the same number of samples in the counts df and the organs metadata df? 
ncol(counts2)==nrow(organs) # TRUE

# Drop the decimals introduced into the column names
#colnames(counts2) <- keeps

# Reorder the columns (samples) in the counts2 df to be in the
# same order as the samples in the organs metadata df (rows)
# Necessary because voom() assumes the rows of the design matrix are the samples
# and they are in the same order as the column of the counts matrix
idx <- match(organs$SAMPID, colnames(counts2))
counts3 <- as.data.table(counts2)
ordered_counts <- counts3[,..idx]
all(organs$SAMPID==colnames(ordered_counts)) # TRUE

# Make design matrix
mat <- model.matrix(~ 0 + organ + SMGEBTCH + SMRIN + SMTSISCH, data=organs)
nrow(mat)==ncol(ordered_counts) # TRUE

#_______________________________________________________________________________ 
# Perform quantile normalization 
#_______________________________________________________________________________ 
# Does the number of columns (samples) in the ordered counts matrix
# match the number of rows (samples) in the organs df?
counts_mat <- as.matrix(ordered_counts)
ncol(counts_mat)==nrow(organs) # TRUE
nrow(ordered_counts) # 56,200

# Apply voom with quality weights  normalization and quantile normalization
#voom_obj <- voom(counts=counts_mat[1:20,1:20], design=mat[1:20,], normalize.method="quantile", save.plot=TRUE) # TEST 
print("line 57")
voom_obj <- voomWithQualityWeights(counts=counts_mat, design=mat, normalize.method="quantile", save.plot=TRUE) # TEST 
print("line 59")

#write.csv(voom_obj$weights, "/scratch/mjpete11/linear_models/data/batch_voom_qnorm_weights.csv")
#gene_weights <- fread("/scratch/mjpete11/linear_models/data/batch_voom_qnorm_weights.csv")

# Perform differential expression
#fit <- lmFit(voom_obj, mat[1:20,])
fit <- lmFit(voom_obj, mat)
print("line 67")
fit <- eBayes(fit)

# Write the logFC and moderated t statistics to file so they can be added to the organs dataframe
write.csv(fit$coefficients, "/scratch/mjpete11/linear_models/data/voomWithArrayWeights_logFC.csv")
write.csv(fit$t, "/scratch/mjpete11/linear_models/data/voomWithArrayWeights_moderatedt.csv")
print("line 73")

# Adding logFC and moderated t statistics will be in a separate script (limma_tstats_organs.R)
