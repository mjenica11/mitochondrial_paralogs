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

# How many samples are there?
ncol(counts2)
length(unique(colnames(counts2))) # 296 samples --> 148 paired heart and liver samples

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
# Perform voom + quantile normalization 
#_______________________________________________________________________________ 
# Does the number of columns (samples) in the ordered counts matrix
# match the number of rows (samples) in the organs df?
counts_mat <- as.matrix(ordered_counts)
ncol(counts_mat)==nrow(organs) # TRUE
nrow(ordered_counts) # 56,200
print("Made ordered count matrix")

# Apply voom with quality weights  normalization and quantile normalization
#voom_obj <- voomWithQualityWeights(counts=counts_mat[1:20,1:20], design=mat[1:20,], normalize.method="quantile") # TEST 
voom_obj <- voom(counts=counts_mat, design=mat, normalize.method="quantile")
print("voom completed")

# Estimate array weights on top of voom weights	
aw <- arrayWeights(voom_obj, design=mat, method="auto", maxiter=maxiter)
print("first weights completed")

# Update voom weights now using the array weights
voom_obj <- voom(counts_mat, design=mat, weights=aw, normalize.method="quantile", plot=FALSE, span=span, ...)
print("second voom completed")

# Update array weights again
aw <- arrayWeights(voom_obj, design=mat, method="auto")
print("second weights completed")

# Incorporate the array weights into the voom weights
voom_obj$weights <- t(aw * t(voom_obj$weights))
voom_obj$targets$sample.weights <- aw

# Write the matrix of normalized expression values on the log2 scale
write.csv(voom_obj$E, "/scratch/mjpete11/linear_models/data/voomWithArrayWeights_matrix1.csv")
print("wrote the E object to file")

# Perform differential expression
#fit <- lmFit(voom_obj, mat[1:20,])
fit <- lmFit(voom_obj, mat)
print("lmFit completed")
fit <- eBayes(fit)
print("eBayes completed")

# Calculate the logFC and confidence intervals
# Set the coefficient to be the contrast between heart and liver
res <- topTable(fit, coef=1, number=Inf, genelist=SLC, adjust.method="BH",
				sort.by="logFC", p.value=0.05, lfc=0, confint=TRUE)
write.csv(res, "/scratch/mjpete11/linear_models/data/voomWithArrayWeights_topTable.1csv")

# Write the logFC and moderated t statistics to file so they can be added to the organs dataframe
write.csv(fit$coefficients, "/scratch/mjpete11/linear_models/data/voomWithArrayWeights_logFC1.csv")
write.csv(fit$t, "/scratch/mjpete11/linear_models/data/voomWithArrayWeights_moderatedt1.csv")
print("Done")

# Kernel density plots are in violin_plots.R
