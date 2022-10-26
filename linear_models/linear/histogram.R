# Plot distribution of SLC25 gene raw counts

# Libraries
library(plyr)
library(limma)
library(data.table)
library(tidyverse)
library(reshape)
library(stringr)
library(ggplot2)
library(ggpubr)

# Read in quantile normalized counts
# Start with 9 samples since I keep getting out of memory errors
counts <- fread("/scratch/mjpete11/linear_models/data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct", sep="\t") 

# List of SLC25 paralogs
SLC <- c("SLC25A1", "SLC25A2", "SLC25A3", "SLC25A4", "SLC25A5", "SLC25A6",
		 "UCP1", "UCP2", "UCP3", "SLC25A10", "SLC25A11", "SLC25A12", 
		 "SLC25A13", "SLC25A14", "SLC25A15", "SLC25A16", "SLC25A17", 
		 "SLC25A18", "SLC25A19", "SLC25A20", "SLC25A21", "SLC25A22", 
		 "SLC25A23", "SLC25A24", "SLC25A25", "SLC25A26", "SLC25A27", 
		 "SLC25A28", "SLC25A29", "SLC25A30", "SLC25A31", "SLC25A32", 
		 "SLC25A33", "SLC25A34", "SLC25A35", "SLC25A36", "SLC25A37", 
		 "SLC25A38", "SLC25A39", "SLC25A40", "SLC25A41", "SLC25A42", 
		 "SLC25A43", "SLC25A44", "SLC25A45", "SLC25A46", "SLC25A47",
		 "MTCH1", "MTCH2", "SLC25A50",  "SLC25A51", "SLC25A52", "SLC25A53")

# Subset the SLC25 genes
counts <- as.data.frame(counts)
sub_df <- counts[counts$"Description" %in% SLC, ]

# Write to file
#write.table(sub_df, "/scratch/mjpete11/linear_models/data/sub_df.csv", sep=",")

# Are any genes missing? 
setdiff(SLC, sub_df$'Description') # No

# Don't think I want to do filtering, since the objective is to look at the
# distributions of the raw counts
# Median filter: Function to filter rows that do not have a median >= condition
#median_filter <- function(DF, thresh){
#		DAT <- DF[which(apply(DF[,-c(1:2)],1,median) > thresh), ]
#		return(DAT)
#}

# Filter rows with a median of < 10 counts
#counts2 <- median_filter(DF=counts, thresh=10)

# Number of genes remaining
#nrow(counts) # 56,200

# Read in sample attributes files
file2 <- read.csv("/scratch/mjpete11/linear_models/data/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt", sep="\t")

# Subset heart left ventricle and liver samples into two separate dfs
heart <- file2[file2$SMTSD %in% "Heart - Left Ventricle", ]
liver <- file2[file2$SMTSD %in% "Liver", ]

# Subset the SLC25 gene count df by the sample IDs that match the IDs in file3 df
heart2 <- sub_df %>% select(contains(heart$SAMPID))
liver2 <- sub_df %>% select(contains(liver$SAMPID))

# Append the gene ensemble ID and common name
heart3 <- cbind('gene'=sub_df$'Description', heart2)
liver3 <- cbind('gene'=sub_df$'Description', liver2)

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
heart4$SMGEBTCH <- file2$SMGEBTCH[match(heart4$SAMPID, file2$SAMPID)]
liver4$SMGEBTCH <- file2$SMGEBTCH[match(liver4$SAMPID, file2$SAMPID)]

heart4$SMGEBTCHT <- file2$SMGEBTCHT[match(heart4$SAMPID, file2$SAMPID)]
liver4$SMGEBTCHT <- file2$SMGEBTCHT[match(liver4$SAMPID, file2$SAMPID)]

heart4$SMRIN <- file2$SMRIN[match(heart4$SAMPID, file2$SAMPID)]
liver4$SMRIN <- file2$SMRIN[match(liver4$SAMPID, file2$SAMPID)]

heart4$SMTSISCH <- file2$SMTSISCH[match(heart4$SAMPID, file2$SAMPID)]
liver4$SMTSISCH <- file2$SMTSISCH[match(liver4$SAMPID, file2$SAMPID)]

# Combine into one df
organs <- rbind(heart4, liver4)

# Add sample ID col to meta
heart4$SUBJID <- str_sub(heart4$SAMPID, start=1L, end=10L)
liver4$SUBJID <- str_sub(liver4$SAMPID, start=1L, end=10L)

# Drop columns 
heart4$SMGEBTCH <- NULL
heart4$SMGEBTCHT <- NULL

liver4$SMGEBTCH <- NULL
liver4$SMGEBTCHT <- NULL

# Dataframe with only samples from people who donated both a heart and liver: organs
# e.g. paired samples
organs_tmp <- merge(heart4, liver4, by="SUBJID") 
length(unique(organs_tmp$SUBJID)) # 148 individuals
tmp1 <- organs_tmp[,c("SUBJID", "gene.x", "organ.x", "SAMPID.x", "value.x", "SMRIN.x", "SMTSISCH.x")]
tmp2 <- organs_tmp[,c("SUBJID", "gene.y", "organ.y", "SAMPID.y", "value.y", "SMRIN.y", "SMTSISCH.y")]
colnames(tmp1) <- c("SUBJID", "gene", "organ", "SUBJID", "value", "SMRIN", "SMTSISCH")
colnames(tmp2) <- c("SUBJID", "gene", "organ", "SUBJID", "value", "SMRIN", "SMTSISCH")
organs <- rbind(tmp1, tmp2)
colnames(organs)[4] <- "SAMPID"

# Write organs df to file
write.table(organs, "/scratch/mjpete11/linear_models/data/raw_counts_organs.csv", sep=",")

# Read organs df back in
#organs <- read.csv("/scratch/mjpete11/linear_models/data/raw_counts_organs.csv", sep=",")

# Historgam function
# Removed the x and y axis because it wouldn't plot half the genes
# even if I set it to the maximum possible value :/
histogram <- function(GENE){
     dat <- organs %>% filter(gene == GENE)
#     dat2 <- distinct(na.omit(dat))
     p <- ggplot(dat, aes(value)) +
       geom_histogram(bins=100, na.rm=TRUE) +
       ylab("frequency") +
       xlab("counts") +
#	   xlim(c(0,320900)) +
#	   ylim(c(0,200)) +
       ggtitle(paste0("Histograms of raw ", GENE, " expression")) 
	 ggsave(paste0("/scratch/mjpete11/linear_models/linear/raw_counts_histograms/", GENE, ".pdf"), p, device="pdf")
}
plts <- Map(histogram, GENE=SLC)


