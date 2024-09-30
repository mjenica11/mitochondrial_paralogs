#!/usr/bin/env Rscript

# Plot umap projections

# Load libraries
library(umap)
library(ggplot2)
library(stringr)
library(edgeR)
library(data.table)
library(tidyverse)
library(readxl)
library(janitor)
library(viridis)

# Read in umap projections 
#proj <- fread("/scratch/mjpete11/mitochondrial_paralogs/sweet_data/umap/umap_normalized.csv", sep=",")
proj <- fread("/scratch/mjpete11/mitochondrial_paralogs/sweet_data/umap/umap_no_qnorm_normalized.csv", sep=",")
#proj <- fread("/scratch/mjpete11/mitochondrial_paralogs/sweet_data/umap/filtered_only_umap_normalized.csv", sep=",")

# Check the dimensions of the objects
dim(proj) # 63 3
#dim(proj) # 63 3
head(proj)

# Read in sample attributes files
manifest <- read.csv("/scratch/mjpete11/mitochondrial_paralogs/sweet_data/count_matrices/metadata.csv", sep=",")

# Read in the second metadata file with the RIN scores
metadata <- read_excel("/scratch/mjpete11/mitochondrial_paralogs/sweet_data/count_matrices/sweet_metadata2.xlsx")

# Drop the 37th sample; this one couldn't be processed 
#class(metadata)
#metadata <- data.frame(metadata)
#metadata <- row_to_names(metadata, row_number=1)
#dim(metadata) # 64 7; 62 8
#head(metadata)
#metadata[c(37),] # DCM46 Female White Not Hispanic/Latino                46                       NA  8.3000000000000007
#metadata <- metadata[-c(37),]
#dim(metadata) # 63 7
#
#metadata[1:5,1:7]
#colnames(metadata)
#metadata$RIN <- as.numeric(metadata[[7]])
#metadata$RIN
#
## Add the entire metadata for the samples in the corresponding umap projections
meta <- manifest[manifest$'BioSample' %in% proj$V1,]
meta[,meta$BioSample=="SAMN09484097"] # This sample is missing; should not be in the projections
meta <- meta[, c(1,2,7,14,25,28)] # Drop irrelevant variables
head(meta)
dim(meta) # 63 6; 64 6

# Check class factor
apply(meta, 2, class)
class(meta)

# Factor response variables
meta$disease <- as.factor(meta$disease)
meta$sex <- as.factor(meta$sex)
meta$AGE <- as.numeric(meta$AGE)

# Check class factor
apply(meta, 2, class)

# Add the RIN
nrow(meta) # 63
nrow(metadata) # 64
#head(metadata$RIN)
#meta$RIN <- metadata$RIN[-c(37)]
#meta$RIN

# Manually check that the labels match the samples
dim(meta) # 63 6
meta[1:10,]
meta[30:40,]

# Write to file and use this for the PCA plots
#write.csv(meta, "/scratch/mjpete11/mitochondrial_paralogs/sweet_data/count_matrices/subset_metadata.csv")

# Rename cols in projection dfs
colnames(proj) <- c("BioSample", "UMAP_1", "UMAP_2")

# Merge the projections with the metadata
proj <- merge(proj, meta, by="BioSample")
head(proj)

#______________________________________________________________________________
# Plots 
#______________________________________________________________________________
# Function to plot color points by one feature
rm(umap_plot)

set.seed(1993) 

umap_plot <- function(PROJ, FEATURE){
#umap_plot <- function(PROJ){
	p = ggplot(PROJ, aes(x=.data[["UMAP_1"]], y=.data[["UMAP_2"]], color = .data[[FEATURE]])) + # aes_string was deprecated and update doesn't work
	#p = ggplot(PROJ, aes(x=.data[["UMAP_1"]], y=.data[["UMAP_2"]], color = .data[["RIN"]])) +
		geom_point(size = 1) +
		theme_classic(base_size = 18) +
		xlab("UMAP 1") +
		ylab("UMAP 2") +
		coord_fixed() +
		theme(axis.ticks=element_blank(), axis.text=element_blank())
	p <- p + labs(color =  FEATURE) 
	#p <- p + labs(color =  "RIN") 
	p <- p + guides(shape = guide_legend(override.aes = list(size = 3))) 
#	p <- p + guides(color = guide_legend(override.aes = list(size = 3))) 
	p <- p + theme(legend.title = element_text(size = 5), 
				   legend.text = element_text(size = 5))
    p <- p + scale_color_viridis() 
	return(p)
}

#umap_plot(PROJ=proj, FEATURE="RIN")
umap_plot(PROJ=proj, FEATURE="disease")

# Subset projections by samples 
#proj1 <- proj[(proj$BioSample %in% meta$BioSample), ]
head(proj)
class(proj)
proj <- data.frame(proj)

# Print plots
# Color projections by disease type 
#ggsave(umap_plot(PROJ = proj, FEATURE = "disease"), file = "/scratch/mjpete11/mitochondrial_paralogs/sweet_data/umap/normalized_disease_umap.pdf")
ggsave(umap_plot(PROJ = proj, FEATURE = "disease"), file = "/scratch/mjpete11/mitochondrial_paralogs/sweet_data/umap/no_qnorm_normalized_disease_umap.pdf")
#ggsave(umap_plot(PROJ = proj, FEATURE = "disease"), file = "/scratch/mjpete11/mitochondrial_paralogs/sweet_data/umap/filtered_disease_umap.pdf")
                                   
# Color projections by RIN 
#ggsave(umap_plot(PROJ = proj, FEATURE = "RIN"), file = "/scratch/mjpete11/mitochondrial_paralogs/sweet_data/umap/no_qnorm_normalized_RIN_umap.pdf")
#ggsave(umap_plot(PROJ = proj, FEATURE = "RIN"), file = "/scratch/mjpete11/mitochondrial_paralogs/sweet_data/umap/filtered_RIN_umap.pdf")
                                   
# Color projections by age 
#ggsave(umap_plot(PROJ = proj, FEATURE = "AGE"), file = "/scratch/mjpete11/mitochondrial_paralogs/sweet_data/umap/no_qnorm_normalized_age_umap.pdf")
ggsave(umap_plot(PROJ = proj, FEATURE = "AGE"), file = "/scratch/mjpete11/mitochondrial_paralogs/sweet_data/umap/no_qnorm_normalized_age_umap.pdf")
#ggsave(umap_plot(PROJ = proj, FEATURE = "AGE"), file = "/scratch/mjpete11/mitochondrial_paralogs/sweet_data/umap/filtered_age_umap.pdf")

# Color projections by sex 
#ggsave(umap_plot(PROJ = proj, FEATURE = "sex"), file = "/scratch/mjpete11/mitochondrial_paralogs/sweet_data/umap/no_qnorm_normalized_sex_umap.pdf")
ggsave(umap_plot(PROJ = proj, FEATURE = "sex"), file = "/scratch/mjpete11/mitochondrial_paralogs/sweet_data/umap/no_qnorm_normalized_sex_umap.pdf")
#ggsave(umap_plot(PROJ = proj, FEATURE = "sex"), file = "/scratch/mjpete11/mitochondrial_paralogs/sweet_data/umap/filtered_sex_umap.pdf")

