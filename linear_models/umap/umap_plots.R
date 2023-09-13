#!/usr/bin/env Rscript

# Plot umap projections

# Load libraries
library(umap)
library(ggplot2)
library(stringr)
library(edgeR)
library(data.table)
library(tidyverse)

# Read in GTEx counts
counts <- fread("/scratch/mjpete11/linear_models/data/combat_batch_adjusted_counts.csv", sep=",")

# Read in sample attributes files
manifest <- read.csv("/scratch/mjpete11/linear_models/data/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt", sep="\t")

# Subset heart left ventricle and liver samples into two separate dfs
heart <- manifest[manifest$SMTSD %in% "Heart - Left Ventricle", ]
liver <- manifest[manifest$SMTSD %in% "Liver", ]

# Subset the SLC25 gene count df by the sample IDs that match the IDs in file3 df
heart2 <- counts %>% select(contains(heart$SAMPID))
liver2 <- counts %>% select(contains(liver$SAMPID))

# Append the gene ensemble ID and common name
heart3 <- cbind('gene'=counts$'Description', heart2)
liver3 <- cbind('gene'=counts$'Description', liver2)

heart3 <- cbind('gene'=counts[,1], heart2)
liver3 <- cbind('gene'=counts[,1], liver2)
colnames(heart3)[1] <- 'gene' 
colnames(liver3)[1] <- 'gene' 

# Add a columns with the organ in each df
heart3$organ <- 'heart' 
liver3$organ <- 'liver' 

# Reshape dataframe so it can be converted to a design matrix object
heart4 <- melt(data = heart3,
			   id.vars = c("gene", "organ"),
			   measure.vars = colnames(heart3)[3:ncol(heart3)-1],
			   variable.name = "samples",
			   value.name = "counts")

liver4 <- melt(data = liver3,
			   id.vars = c("gene", "organ"),
			   measure.vars = colnames(liver3)[3:ncol(liver3)-1],
			   variable.name = "samples",
			   value.name = "counts")

# Change the name of the 'variable' column to 'SAMPID' to match columns
colnames(heart4)[3] <- "SAMPID" 
colnames(liver4)[3] <- "SAMPID" 

# Add column with the expression batch ID (SMGEBTCH) and the type of genotype
# or expression batch (SMGEBTCHT)
heart4$SMGEBTCH <- manifest$SMGEBTCH[match(heart4$SAMPID, manifest$SAMPID)]
liver4$SMGEBTCH <- manifest$SMGEBTCH[match(liver4$SAMPID, manifest$SAMPID)]

heart4$SMGEBTCHT <- manifest$SMGEBTCHT[match(heart4$SAMPID, manifest$SAMPID)]
liver4$SMGEBTCHT <- manifest$SMGEBTCHT[match(liver4$SAMPID, manifest$SAMPID)]

heart4$SMRIN <- manifest$SMRIN[match(heart4$SAMPID, manifest$SAMPID)]
liver4$SMRIN <- manifest$SMRIN[match(liver4$SAMPID, manifest$SAMPID)]

heart4$SMTSISCH <- manifest$SMTSISCH[match(heart4$SAMPID, manifest$SAMPID)]
liver4$SMTSISCH <- manifest$SMTSISCH[match(liver4$SAMPID, manifest$SAMPID)]

# Combine into one df
organs <- rbind(heart4, liver4)

# Factor response variable
organs$organ <- as.factor(organs$organ)

#______________________________________________________________________________
# Plots 
#______________________________________________________________________________
# Read in umap projections
proj <- read.csv("/scratch/mjpete11/linear_models/results2/umap_combat_normalized.csv", sep=",", header=TRUE)

# Rename cols in projection dfs
colnames(proj) <- c("SAMPID", "UMAP_1", "UMAP_2")

# Add metadata to projections for plotting; bind using sample ID
proj2 <- merge(proj, manifest, by = "SAMPID")

# Function to plot color points by one feature
umap_plot <- function(PROJ, FEATURE){
	p = ggplot(PROJ, aes_string("UMAP_1", "UMAP_2", color = FEATURE)) +
		geom_point(size = 0.5) +
		theme_classic(base_size = 18) +
		xlab("UMAP 1") +
		ylab("UMAP 2") +
		coord_fixed() +
		theme(axis.ticks=element_blank(), axis.text=element_blank())
	p <- p + labs(fill =  FEATURE) 
	p <- p + guides(shape = guide_legend(override.aes = list(size = 3))) 
	p <- p + guides(color = guide_legend(override.aes = list(size = 3))) 
	p <- p + theme(legend.title = element_text(size = 0), 
				   legend.text = element_text(size = 0))
	return(p)
}

# Subset projections by heart and liver 
proj3 <- proj2[(proj2$SAMPID %in% organs$SAMPID), ]

# Drop rows if NA in columns of interest (sample ID, tissue type, ischemic time, batch ID, and RIN number)
proj3 <- proj3 %>% drop_na(SAMPID, SMTSD, SMTSISCH, SMGEBTCH, SMRIN) 

# Print plots
# Color projections by tissue type 
ggsave(umap_plot(PROJ = proj3, FEATURE = "SMTSD"), file = "/scratch/mjpete11/linear_models/results2/combat_normalized_tissue_umap.pdf")
                                   
# Color projections by batch 
ggsave(umap_plot(PROJ = proj3, FEATURE = "SMGEBTCH"), file = "/scratch/mjpete11/linear_models/results2/combat_normalized_batch_umap.pdf")
                                   
# Color projections by ischemic time 
ggsave(umap_plot(PROJ = proj3, FEATURE = "SMTSISCH"), file = "/scratch/mjpete11/linear_models/results2/combat_normalized_ischemic_time_umap.pdf")

# Color projections by RIN 
ggsave(umap_plot(PROJ = proj3, FEATURE = "SMRIN"), file = "/scratch/mjpete11/linear_models/results2/combat_normalized_RIN_umap.pdf")

