# Purpose: Linear models to asses if SLC paralogs are dependent on energetic state?

# Libraries
library(plyr)
library(limma)
library(edgeR)
library(data.table)
library(tidyverse)
library(reshape)
library(stringr)
library(moderndive)
library(infer)
library(performance)
library(see)
library(ggplot2)
library(ggpubr)

# Read in GTEx counts
counts <- fread("/scratch/mjpete11/linear_models/data/quantile_normalized_counts.csv", sep=",")

# Add gene name column
gene_names <- fread("/scratch/mjpete11/linear_models/data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct", sep="\t",
					select=c("Description"))
counts$"Description" <- gene_names 

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
		 "SLC25A48", "MTCH1", "MTCH2", "SLC25A51", "SLC25A52", "SLC25A53")

# Subset the SLC25 genes
sub_df <- counts[counts$"Description" %in% SLC, ]

# Only 49 genes are present; which are missing?
setdiff(SLC, sub_df$'Description') 
# "SLC25A7"  "SLC25A8"  "SLC25A9"  "SLC25A49" "SLC25A50"
# I searched for the ENSG ID but I couldn't find them
# Update: tried the following aliases found in the GTEx portal: UCP1, UCP2,
# UCP3, MTCH1, MTCH2

# Median filter: Function to filter rows that do not have a median >= condition
median_filter <- function(DF, thresh){
		DAT <- DF[which(apply(DF[,-c(1:2)],1,median) > thresh), ]
		return(DAT)
}

# Filter rows with a median of < 10 counts
#counts <- median_filter(DF=counts, thresh=10)
counts <- median_filter(DF=counts, thresh=10)

# Number of genes remaining
nrow(counts) # 19,514

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
#organs <- rbind(heart4, liver4)
#organs <- left_join(heart4, liver4, by="gene")

# Add a column with the individual IDs
# Add sex and age to manifest
# Read in separate file
meta <- read.csv("/scratch/mjpete11/linear_models/data/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt", header=TRUE, sep = "\t")

# Add sample ID col to meta
heart4$SUBJID <- str_sub(heart4$SAMPID, start=1L, end=10L)
liver4$SUBJID <- str_sub(liver4$SAMPID, start=1L, end=10L)

# Subset to only heart and liver samples that are from the same person
heart4$SMGEBTCH <- NULL
#heart4$SMRIN <- NULL
#heart4$SMTSISCH <- NULL
heart4$SMGEBTCHT <- NULL

liver4$SMGEBTCH <- NULL
#liver4$SMRIN <- NULL
#liver4$SMTSISCH <- NULL
liver4$SMGEBTCHT <- NULL

# Dataframe with only samples from people who donated both a heart and liver: organs
# e.g. paired samples
organs_tmp <- merge(heart4, liver4, by="SUBJID") # 148 samples
tmp1 <- organs_tmp[,c("SUBJID", "gene.x", "organ.x", "SAMPID.x", "value.x", "SMRIN.x", "SMTSISCH.x")]
tmp2 <- organs_tmp[,c("SUBJID", "gene.y", "organ.y", "SAMPID.y", "value.y", "SMRIN.y", "SMTSISCH.y")]
colnames(tmp1) <- c("SUBJID", "gene", "organ", "SAMPID", "value", "SMRIN", "SMTSISCH")
colnames(tmp2) <- c("SUBJID", "gene", "organ", "SAMPID", "value", "SMRIN", "SMTSISCH")
organs <- rbind(tmp1, tmp2)

# Write organs df to file
#write.csv(organs, "/scratch/mjpete11/linear_models/linear/organs.csv", sep=",")

# Read organs df back in
organs <- read.csv("/scratch/mjpete11/linear_models/linear/organs.csv", sep=",")

# Drop the index columns
organs$X <- NULL

# Regress out effect of RIN and ischemic time and remove from the model 
# before making the violin plots

# Subset one gene from the organ df
A1 <- organs %>% filter(gene == "SLC25A1")
dA1 <- distinct(na.omit(A1))
# Double check that for each value in the subject ID column there is a 
# corresponding heart and liver value in the organ column
# e.g. How many unique values in organ for each value of subject ID after
# subsetting by gene? Are any unique values not equal to 2?
any(tapply(dA1$organ, dA1$SUBJID, function(x) length(unique(x))) != 2) #FALSE
# Regress out the effect due to RIN and ischemic time
model1 <- lm(value ~ organ + SMRIN + SMTSISCH, data=dA1)
# Extract residuals from model and add to dataframe
dA1$corrected <- residuals(model1)
# Are the number of rows divided by 2 equal to the number of unique subject IDs?
nrow(dA1)/2==148 # TRUE

# Test violin jitter plot function
p <- ggplot(dA1, aes(x = organ, y = corrected)) +
	 geom_violin(aes(colour = organ)) +
	 geom_jitter(aes(colour = organ))  +
	 ylab("value") +
	 xlab("organ") +
	 ggtitle(paste0("Jitter plot of SLC25A1 expression between heart and liver")) +
	 stat_compare_means(method = "t.test") +
	 stat_summary(fun.data = "mean_cl_boot", geom = "crossbar", colour = "red", width = 0.2)
ggsave(paste0("SLC25A1_violin_plot.pdf"))

# Violin and jitter plot function
violin_plots <- function(GENE){
			dat <- organs %>% filter(gene == GENE)
			dat2 <- na.omit(dat)
			dat3 <- distinct(na.omit(dat2))
			any(tapply(dat3$organ, dat3$SUBJID, function(x) length(unique(x))) != 2) #FALSE
			model1 <- lm(value ~ organ + SMRIN + SMTSISCH, data=dat3)
			dat3$corrected <- residuals(model1)
			nrow(dat3)/2==148 # TRUE
			p <- ggplot(dat3, aes(x = organ, y = value)) +
			geom_violin(aes(colour = organ)) +
			geom_jitter(aes(colour = organ))  +
			ylab("value") +
			xlab("organ") +
			ggtitle(paste0("Violin plot of ", GENE, " expression between heart and liver")) +
			stat_compare_means(method="t.test") +
			stat_summary(fun.data = "mean_cl_boot", geom = "crossbar", colour = "red", width = 0.2)
			ggsave(paste0(GENE, "_violin_plot.pdf"))
}
plots <- Map(violin_plots, GENE=SLC)
