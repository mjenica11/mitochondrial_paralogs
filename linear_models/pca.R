# PCA plots of unadjusted, quantile normalize, and combat normalized GTEx count data

library(data.table)
library(stats)
library(ggplot2)
library(ggfortify)
library(tidyverse)

# Read in quantile normalized GTEx data
counts <- fread("/scratch/mjpete11/linear_models/data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct", sep="\t")

# Read in metadata with batch variables
manifest <- read.csv("/scratch/mjpete11/linear_models/data/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt", header=TRUE, sep = "\t")

# Remove all rows with NA in either columns
manifest <- manifest %>% drop_na(SAMPID, SMTSD, SMTSISCH, SMGEBTCH, SMRIN) 

# Median filter: Function to filter rows that do not have a median >= condition
median_filter <- function(DF, thresh){
		DF <- DF[which(apply(DF[,-c(1:2)],1,median) > thresh), ]
		return(DF)
}

# Filter rows with a median of < 10 counts
counts <- median_filter(DF=counts, thresh=10)

# Drop samples from manifest that are missing in the count data
#rownames(manifest) <- manifest$SAMPID
manifest <- manifest[(manifest$SAMPID %in% colnames(counts)[3:ncol(counts)]),]
samples <- manifest$SAMPID

# Only keep samples that are present in the manifest 
counts <- counts[,..samples]

# PCA on unnormalized data
#pca_obj <- prcomp(counts[,3:ncol(counts)])
pca_obj <- prcomp(counts[,3:20])

# convert to dataframe
pca_df <- as.data.frame(pca_obj[2]$rotation)

# Combine pca df with metadata for plotting
pca_df$SAMPID <- rownames(pca_df)  
pca_df <- merge(pca_df, manifest, by="SAMPID")

#as above, create a PCA plot for comparison to the uncorrected data
#cols <- c("UHR" = "#481567FF", "HBR" = "#1F968BFF")
#p1 = ggplot(data=pca_df, aes(x=PC1, y=PC2))
#p1 = p1 + geom_point(size=3)
#p1 = p1 + stat_ellipse(type="norm", linetype=2)
#p1 = p1 + labs(title="PCA of unnormalized GTEx counts")
#p1 = p1 + scale_colour_manual(values = cols)

pdf(file="/scratch/mjpete11/linear_models/results2/unnormalized_PCA_2.pdf")
autoplot(pca_obj, data=pca_df, colour="SMGEBTCH")
dev.off()
