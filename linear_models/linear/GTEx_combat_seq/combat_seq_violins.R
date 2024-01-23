# Purpose: Linear models to asses if SLC paralogs are dependent on energetic state?

# Libraries
library(plyr)
library(limma)
library(data.table)
library(tidyverse)
library(reshape)
library(stringr)
library(ggplot2)
library(ggpubr)
library(edgeR)

# Read in quantile normalized counts
counts <- fread("/scratch/mjpete11/linear_models/data/combat_seq_filtered.csv", sep=",") # float

# Drop the index column
counts$V1 <- NULL

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
counts <- as.data.frame(counts)
sub_df <- counts[counts$"Description" %in% SLC, ]

# Only 49 genes are present; which are missing?
setdiff(SLC, sub_df$'Description') 

# Number of genes remaining
nrow(counts) # 51,259 

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

# Drop the extra row 
heart3 <- heart3[-c(54),]
liver3 <- liver3[-c(54),]

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

# Add sample ID col to meta
heart4$SUBJID <- str_sub(heart4$SAMPID, start=1L, end=10L)
liver4$SUBJID <- str_sub(liver4$SAMPID, start=1L, end=10L)

# Dataframe with only samples from people who donated both a heart and liver: organs
# e.g. paired samples
organs_tmp <- merge(heart4, liver4, by="SUBJID") 
length(unique(organs_tmp$SUBJID)) # 148 individuals

# Split wide df into separate dfs so they can be recombined into a long format df
tmp1 <- organs_tmp[,c("SUBJID", "gene.x", "organ.x", "SAMPID.x", "value.x")]
tmp2 <- organs_tmp[,c("SUBJID", "gene.y", "organ.y", "SAMPID.y", "value.y")]

# Rename columns to facilitate binding
colnames(tmp1) <- c("SUBJID", "gene", "organ", "SAMPID", "value")
colnames(tmp2) <- c("SUBJID", "gene", "organ", "SAMPID", "value")

# Are there any duplicated rows?
any(duplicated(tmp1)) # TRUE
any(duplicated(tmp2)) # TRUE

# Drop duplicate rows
tmp1_2 <- tmp1[!duplicated(tmp1),]
tmp2_2 <- tmp2[!duplicated(tmp2),]

# Are there any duplicated rows?
any(duplicated(tmp1_2)) # FALSE 
any(duplicated(tmp2_2)) # FALSE

# Add row number column to facilitate vertical binding
tmp1_2$index <- seq(from=1, to=nrow(tmp1_2), by=1) 
tmp2_2$index <- seq(from=1, to=nrow(tmp2_2), by=1) 

# Do the row indices exactly match?
all(tmp1_2$index==tmp2_2$index)==TRUE # TRUE

# Do the heart and liver dfs have the same number of rows?
nrow(tmp1_2)==nrow(tmp2_2) # TRUE

# Delete some old object to open space
rm(counts)
rm(tmp1)
rm(tmp2)

# Bind liver and heart dfs vertically 
organs <- rbind(tmp1_2, tmp2_2)

#organs <- merge(tmp1, tmp2, by=c("SUBJID"))
nrow(organs)==nrow(tmp1_2)+nrow(tmp2_2) # TRUE
any(duplicated(organs)) # FALSE

# Expected number of rows?
148 * 53 * 2 # 15,688
nrow(organs)==15688 # TRUE
any(duplicated(organs)) # FALSE

# Convert the counts to log2(CPM) scale
organs$log2_cpm <- cpm(organs$value, log=TRUE, prior.count=0.5)

# Write organs df to file (02-20-2022)
#write.table(organs, "/scratch/mjpete11/linear_models/data/organs.csv", sep=",")

# Read organs df back in
#organs <- read.csv("/scratch/mjpete11/linear_models/data/organs.csv", sep=",")

# Remove samples >6 standard deviations away
median(organs$log2_cpm) # 4.00
mean(organs$log2_cpm) # 3.32
sd(organs$log2_cpm) # 3.67
sd(organs$log2_cpm) * 6 # +/- 22.1 
above <- 3.43 + 20.6 # 24.03
below <- 3.43 - 20.6 # -17.17

# Are there any samples outside of this range?
range(organs$log2_cpm) # -6.72 to 12.1 

# Number of samples outside of this range
outlier_above <- organs[organs$log2_cpm > 24.03,] # 0 sample  
outlier_below <- organs[organs$log2_cpm < -17.77,] # 0 samples 
nrow(outlier_above);nrow(outlier_below)

# Details for plots
range(organs$log2_cpm)
rm(violin)
rm(plots)

# Function to plot violin plots
violin <- function(GENE){
		 dat <- organs %>% filter(gene==GENE)
	     p <- ggplot(dat, aes(x = organ, y = log2_cpm, fill = organ)) +
					 stat_compare_means(method = "wilcox.test", 
										aes(label = paste("adj.p_value =", after_stat(!!str2lang("p.adj"))*53)), 
										label.x = 1.25, 
										label.y = max(dat[["log2_cpm"]]) + 1.5,
										paired = TRUE) +
					geom_violin(trim = FALSE) +
					stat_summary(fun.data = "mean_sdl", geom="crossbar", width=0.2, alpha=0.1) +
					scale_fill_manual(values = c("lightgreen", "purple")) +
					geom_jitter(size = 1, alpha = 0.9) +
	   		  		labs(x = "organ", y = "log2(CPM(prior.count=0.5))", fill = "") +
	   		  		ggtitle(paste0("Violin plot of ", GENE, " expression with combat-seq correction")) 
	   		  ggsave(paste0("/scratch/mjpete11/linear_models/linear/combat_seq_violins/plots/", GENE, ".png"), device="png")
}
plots <- Map(violin, GENE=SLC)

