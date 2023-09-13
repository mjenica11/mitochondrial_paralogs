#!/bin/bash Rscript

# Perform quantile normalization on count data

# Libraries
library(data.table)
library(tidyverse)

#_______________________________________________________________________________ 
# Make matrix of genes, samples, counts, organs, and batch effects
#_______________________________________________________________________________ 
# Read in GTEx manifest
counts <- fread("/scratch/mjpete11/linear_models/data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct", sep = "\t")

# Read in GTEx counts
normalized <- fread("/scratch/mjpete11/linear_models/data/quantile_normalized_counts.csv", sep=",")

#_______________________________________________________________________________ 
# Density plots
#_______________________________________________________________________________ 
# Density plot of unquantile normalized heart
edata = log2(counts[,3:ncol(counts)] + 1)
edata = edata[rowMeans(edata) > 3, ]
edata = as.matrix(edata)
colramp = colorRampPalette(c(3,"white",2))(20)
pdf("/scratch/mjpete11/linear_models/results2/not_quantile_adjusted_density_plot.pdf")
plot(density(edata[,1]),col=colramp[1],lwd=3,ylim=c(0,1))
for(i in 2:20){lines(density(edata[,i]),lwd=3,col=colramp[i])}
dev.off()

# Density plot of quantile normalized heart
adata = log2(normalized + 1)
adata = adata[rowMeans(adata) > 3, ]
adata = as.matrix(adata)
colramp = colorRampPalette(c(3,"white",2))(20)
pdf("/scratch/mjpete11/linear_models/results2/quantile_adjusted_density_plot.pdf")
plot(density(adata[,1]),col=colramp[1],lwd=3,ylim=c(0,1))
for(i in 2:20){lines(density(adata[,i]),lwd=3,col=colramp[i])}
dev.off()
