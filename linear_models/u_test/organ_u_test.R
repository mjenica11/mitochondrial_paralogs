# Purpose: Linear models to asses if SLC paralogs are dependent on energetic state?

# Libraries
library("ggpubr")
library(data.table)
library(tidyverse)
library(reshape)
library(dplyr)
library(gridExtra)

# Read in GTEx manifest
manifest <- read.csv("/scratch/mjpete11/linear_models/data/sample.tsv", header=TRUE, sep = "\t")

# Make dataframe with sample id, tissue type
# All of the rin number and ischemic time values were missing...
df1 <- data.frame(manifest$"dbgap_sample_id", manifest$"tissue_type")

# Remove all rows with NA in either columns
df2 <- df1[complete.cases(df1), ]

# Read in GTEx counts
counts <- fread("/scratch/mjpete11/linear_models/data/combat_batch_adjusted_counts.csv", sep=",")
unnorm <- fread("/scratch/mjpete11/linear_models/data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct", sep="\t")

# Add gene names to beginning of count df
counts$'Description' <- unnorm$'Description'


# List of SLC25 paralogs
SLC <- c("SLC25A1", "SLC25A2", "SLC25A3", "SLC25A4", "SLC25A5", "SLC25A6", 
		 "SLC25A7", "SLC25A8", "SLC25A9", "SLC25A10", "SLC25A11", "SLC25A12", 
		 "SLC25A13", "SLC25A14", "SLC25A15", "SLC25A16", "SLC25A17", 
		 "SLC25A18", "SLC25A19", "SLC25A20", "SLC25A21", "SLC25A22", 
		 "SLC25A23", "SLC25A24", "SLC25A25", "SLC25A26", "SLC25A27", 
		 "SLC25A28", "SLC25A29", "SLC25A30", "SLC25A31", "SLC25A32", 
		 "SLC25A33", "SLC25A34", "SLC25A35", "SLC25A36", "SLC25A37", 
		 "SLC25A38", "SLC25A39", "SLC25A40", "SLC25A41", "SLC25A42", 
		 "SLC25A43", "SLC25A44", "SLC25A45", "SLC25A46", "SLC25A47",
		 "SLC25A48", "SLC25A49", "SLC25A50", "SLC25A51", "SLC25A52", "SLC25A53")

# Subset the SLC25 genes
sub_df <- counts[counts$"Description" %in% SLC, ]

# Only 49 genes are present; which are missing?
setdiff(SLC, sub_df$'Description') 
# "SLC25A7"  "SLC25A8"  "SLC25A9"  "SLC25A49" "SLC25A50"
# I searched for the ENSG ID but I couldn't find them

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

# Combine into one df
organs <- rbind(heart4, liver4)
#organs <- left_join(heart4, liver4, by="gene")

# Factor response variable
#organs$organ <- as.factor(organs$organ)

# Function to make box plots of each SLC ~ organ individually
# List of subset dfs
SLC_dfs <- c("A1", "A2", "A3", "A4", "A5", "A6", "A7", "A10",
			 "A11", "A12", "A13", "A14", "A15", "A16", "A17", "A18", "A19",
			 "A20", "A21", "A22", "A23", "A24", "A25", "A26", "A27", "A28",
			 "A29", "A30", "A31", "A32", "A33", "A34", "A35", "A36", "A37",
			 "A38", "A39", "A40", "A41", "A42", "A43", "A44", "A45", "A46",
			 "A47", "A48", "A51", "A52", "A53")

SLCs <- c("SLC25A1", "SLC25A2", "SLC25A3", "SLC25A4", "SLC25A5", "SLC25A6", 
		 "SLC25A7", "SLC25A10", "SLC25A11", "SLC25A12", 
		 "SLC25A13", "SLC25A14", "SLC25A15", "SLC25A16", "SLC25A17", 
		 "SLC25A18", "SLC25A19", "SLC25A20", "SLC25A21", "SLC25A22", 
		 "SLC25A23", "SLC25A24", "SLC25A25", "SLC25A26", "SLC25A27", 
		 "SLC25A28", "SLC25A29", "SLC25A30", "SLC25A31", "SLC25A32", 
		 "SLC25A33", "SLC25A34", "SLC25A35", "SLC25A36", "SLC25A37", 
		 "SLC25A38", "SLC25A39", "SLC25A40", "SLC25A41", "SLC25A42", 
		 "SLC25A43", "SLC25A44", "SLC25A45", "SLC25A46", "SLC25A47",
		 "SLC25A48", "SLC25A51", "SLC25A52", "SLC25A53")

models <- sprintf("model%d",1:49)
results <- sprintf("res%d",1:49)

# Boxplot of all SLC genes in heart and liver
group_by(organs, organ) %>%
		summarise(
				  count = n(),
				  median = median(value, na.rm = TRUE),
				  IQR = IQR(value, na.rm = TRUE))

pdf("boxplot.pdf")
p <- ggboxplot(organs, x = "organ", y = "value",
		 	   color = "organ", palette = c("#FFA500", "#FF0000"),
		 	   ylab = "counts", xlab = "organ")
p + stat_compare_means()
dev.off()

# Subset one gene from the organ df
A1 <- organs %>% filter(gene == "SLC25A1")
A1 <- na.omit(A1)
group_by(A1, organ) %>%
		summarise(
				  count = n(),
				  median = median(value, na.rm = TRUE),
				  IQR = IQR(value, na.rm = TRUE))
pdf("A1_boxplot.pdf")
p <- ggboxplot(A1, x = "organ", y = "value",
		 	   color = "organ", palette = c("#FFA500", "#FF0000"),
		 	   ylab = "counts", xlab = "organ")
p + stat_compare_means()
dev.off()

# Function
SLCs <- as.list(SLCs)
box_plots <- function(GENE){
	dat <- organs %>% filter(gene == GENE)
	dat <- na.omit(dat)
	p <- ggboxplot(dat, x = "organ", y = "value",
		 	   color = "organ", palette = c("#FFA500", "#FF0000"),
		 	   ylab = "counts", xlab = "organ")
	p + stat_compare_means()
	ggsave(paste0(GENE, "_boxplot.pdf"))
}
plots <- Map(box_plots, GENE=SLCs)

