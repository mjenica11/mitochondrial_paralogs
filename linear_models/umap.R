#!/usr/bin/env Rscript

# Script to plot the umap projections

# Load libraries
library(umap)
library(ggplot2)
library(stringr)
library(edgeR)
library(data.table)

# Read in counts
counts <- fread("/scratch/mjpete11/linear_models/data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct", sep = "\t")

# Read in annotation file
#meta <- read.table("/scratch/mjpete11/linear_models/data/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt", header = TRUE, sep = "\t")

# Read in file with sex and age
#meta2 <- read.csv("/scratch/mjpete11/linear_models/data/2022.05_FullPhenotypeFile.csv", header = TRUE, sep = ",")

#_______________________________________________________________________________
# umap
#_______________________________________________________________________________
# Convert gene_counts to matrix and transpose for use w/ umap
counts <- counts[,3:ncol(counts)]
res_matrix <- as.matrix(t(counts))

# Apply umap to an expression matrix: genes as rows and samples as columns.
# processed count matrix with sex chr
start_t <- Sys.time()
res_umap <- umap(res_matrix, n_neighbors = 50, min_dist = 0.5)
end_t <- Sys.time()
print(paste("Elapsed time to do umap", end_t - start_t))

# Make df of projections from umap object
res_proj <- data.frame(res_umap$layout)

#_______________________________________________________________________________
# Sanity check 
#_______________________________________________________________________________
# Number of samples in meta
#nrow(meta) # 2,146

# Number of samples in projection
#nrow(res_proj) == nrow(meta) # TRUE

# Check that the samples in the projections are in the same order as in meta
#identical(row.names(res_proj), meta$Sample_ID) # TRUE

# Write umap objects to file
write.csv(res_proj, file = "/scratch/mjpete11/linear_models/results2/umap_unnormalized.csv", row.names = TRUE)
