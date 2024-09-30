#!/usr/bin/env Rscript

# Script to plot the umap projections

# Load libraries
library(umap)
library(edgeR)
library(data.table)
library(dplyr)
library(purrr)
library(janitor)

# Read in counts
#counts <- fread("/scratch/mjpete11/mitochondrial_paralogs/sweet_data/qnorm_voom_normalize/voom_qnorm_counts.csv", sep = ",")
counts <- fread("/scratch/mjpete11/mitochondrial_paralogs/sweet_data/qnorm_voom_normalize/voom_only_counts.csv", sep = ",")
dim(counts) # 24209    65

# Read in filtered not not normalized data 
# Read in clinical data 
#nf_counts <- read.table("/scratch/mjpete11/mitochondrial_paralogs/sweet_data/count_matrices/filtered_non_failing_controls.csv", sep=",")
#icm_counts <- read.table("/scratch/mjpete11/mitochondrial_paralogs/sweet_data/count_matrices/filtered_ischemic_cardiomyopathy.csv", sep=",")
#dcm_counts <- read.table("/scratch/mjpete11/mitochondrial_paralogs/sweet_data/count_matrices/filtered_dilated_cardiomyopathy.csv", sep=",")
#
## Combine into one dataframe
#head(nf_counts) 
#head(icm_counts) 
#head(dcm_counts) 
#
## Dimensions
#dim(nf_counts) # 24148 16 --> 14
#dim(icm_counts) # 25100 15 --> 15 
#dim(dcm_counts) # 25100 38 --> 36 
#
## Combine into one dataframe
#counts <- list(nf_counts, icm_counts, dcm_counts) %>% reduce(inner_join, by="V1")
counts[1:5,1:5]
#head(counts)
#dim(counts) # 24148 67 --> 64 samples 
#colnames(counts)

# Make the first row into the colnames
#counts <- counts %>% row_to_names(row_number=1) 

# Drop the extra 2 gene name columns
#counts <- counts[,-c(17,31)]
#colnames(counts)
dim(counts) # 24148 65; 63 samples + 2 ID columns 

# One of the samples should be missing bc Salmon couldn't process it...
counts$'SAMN09484097' # This sample should be missing and it is 

# Subset to the top 500 most highly expressed genes
#dim(counts) # 24147 66
counts[1:5, 1:5]
class(as.data.frame(counts))
apply(counts, 2, class)
counts_numeric <- apply(counts[,3:ncol(counts)], 2, as.numeric)
apply(counts_numeric, 2, class)
class(counts_numeric)
counts_numeric[1:5,1:5]
min(colMeans(counts_numeric)) # 828
max(colMeans(counts_numeric)) # 1464 
#select = order(rowMeans(counts_numeric[,3:ncol(counts_numeric)]), decreasing=TRUE)[1:500]
#highexprgenes_counts <- counts_numeric[,3:ncol(counts_numeric)][select,]
select = order(rowMeans(counts_numeric), decreasing=TRUE)[1:500]
highexprgenes_counts <- counts_numeric[select,]

#dim(highexprgenes_counts) # 500 64 
dim(highexprgenes_counts) # 500 63 
apply(highexprgenes_counts, 2, class)
highexprgenes_counts[1:5, 1:5]

# Convert gene_highexprgenes_counts to matrix and transpose for use w/ umap
#highexprgenes_counts <- highexprgenes_counts[,3:ncol(highexprgenes_counts)]
res_matrix <- as.matrix(t(highexprgenes_counts))
head(res_matrix)
any(is.na(res_matrix))==TRUE
dim(res_matrix) # 62 500
# Drop the filtered columns
#res_matrix <- na.omit(res_matrix)
#any(is.na(res_matrix))==FALSE

# Set the seed so the results are reproducible
set.seed(1993)

# Apply umap to an expression matrix: genes as rows and samples as columns.
# processed count matrix with sex chr
start_t <- Sys.time()
res_umap <- umap(res_matrix, metric="manhattan", n_neighbors = 3, min_dist = 0.5)
end_t <- Sys.time()
print(paste("Elapsed time to do umap", end_t - start_t))

# Make df of projections from umap object
res_proj <- data.frame(res_umap$layout)

# Check if the dimensions match
dim(res_proj) # filtered data: 63 2; columns are the sample names
head(res_proj)

# Write umap objects to file
#write.csv(res_proj, file = "/scratch/mjpete11/mitochondrial_paralogs/sweet_data/umap/umap_normalized.csv", row.names = TRUE)
write.csv(res_proj, file = "/scratch/mjpete11/mitochondrial_paralogs/sweet_data/umap/umap_no_qnorm_normalized.csv", row.names = TRUE)
#write.csv(res_proj, file = "/scratch/mjpete11/mitochondrial_paralogs/sweet_data/umap/filtered_only_umap_normalized.csv", row.names = TRUE)

