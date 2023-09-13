#!/usr/bin/env Rscript

# Script to plot the umap projections

# Load libraries
library(umap)
library(edgeR)
library(data.table)

# Read in counts
counts <- fread("/scratch/mjpete11/linear_models/data/combat_batch_adjusted_counts.csv", sep = ",")

# Convert gene_counts to matrix and transpose for use w/ umap
counts <- counts[,2:ncol(counts)]
res_matrix <- as.matrix(t(counts))

# Set the seed so the results are reproducible
set.seed(1993)

# Apply umap to an expression matrix: genes as rows and samples as columns.
# processed count matrix with sex chr
start_t <- Sys.time()
res_umap <- umap(res_matrix, n_neighbors = 50, min_dist = 0.5)
end_t <- Sys.time()
print(paste("Elapsed time to do umap", end_t - start_t))

# Make df of projections from umap object
res_proj <- data.frame(res_umap$layout)

# Write umap objects to file
write.csv(res_proj, file = "/scratch/mjpete11/linear_models/results2/umap_combat_normalized.csv", row.names = TRUE)
