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

# Read in voom + qnorm adjusted filtered GTEx count matrix
#voom_matrix <- fread("/scratch/mjpete11/linear_models/data/simulated_data_voom_Eobject.csv", sep=",") 
voom_matrix <- fread("/scratch/mjpete11/linear_models/data/batch_voom_qnorm_matrix.csv", sep=",") 
voom_matrix[1:10,1:10]
dim(voom_matrix) # 51259 15984

# Read in sample attributes files
file2 <- read.csv("/scratch/mjpete11/linear_models/data/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt", sep="\t")

# Drop the index column
voom_matrix$V1 <- NULL

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

# Convert matrix to data.frame
voom_df <- as.data.frame(voom_matrix)
# Add gene name column
gene_names <- fread("/scratch/mjpete11/linear_models/data/filtered_counts.csv", sep=",",
					select=c("Description"))
# Convert gene_names from df to vector
gene_names <- gene_names[["Description"]] 
class(gene_names)
dim(voom_df) # 51259 15984
length(gene_names) # 51259

# Add the gene names (hugo_IDs) column to the voom adjusted counts dataframe
voom_df$"hugo_IDs" <- gene_names 

# Subset the SLC25 genes
sub_df <- voom_df[voom_df$"hugo_IDs" %in% SLC, ]
sub_df[1:5,1:5]

# Number of genes remaining
dim(voom_df) # 56200 15985
dim(sub_df) # 53 15985 

# Keep only unique colnames
voom_df <- subset(sub_df, select=unique(colnames(voom_df)))
dim(voom_df) # 53 297 (148 samples * 2 conditions + 1 gene name columns)

# Move the gene names column to the front
voom_df <- voom_df %>% select("hugo_IDs", everything())
voom_df[1:5,1:5]

# Subset heart left ventricle and liver samples into two separate dfs
heart <- file2[file2$SMTSD %in% "Heart - Left Ventricle", ]
liver <- file2[file2$SMTSD %in% "Liver", ]

# Subset the voom + qnorm adjusted matrix by the sample IDs that match the IDs in file3 df
heart2 <- voom_df %>% select(contains(heart$SAMPID))
liver2 <- voom_df %>% select(contains(liver$SAMPID))

# How many heart/liver samples are there?
ncol(heart2) # 148
ncol(liver2) # 148

# Append the gene ensemble ID and common name
heart3 <- cbind('gene'=voom_df$'hugo_IDs', heart2)
liver3 <- cbind('gene'=voom_df$'hugo_IDs', liver2)

# Add a columns with the organ in each df
heart3$organ <- 'heart' 
liver3$organ <- 'liver' 

# Drop the extra row
heart3 <- heart3[-c(54),]
liver3 <- liver3[-c(54),]

# Reshape dataframe so that gene, organ, sample ID and voom adjusted values
# are columns
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

# Add sample ID col to the long-format datatframes
heart4$SUBJID <- str_sub(heart4$SAMPID, start=1L, end=10L)
liver4$SUBJID <- str_sub(liver4$SAMPID, start=1L, end=10L)

# Dataframe with only samples from people who donated both a heart and liver: organs
# e.g. paired samples
plotting_df_tmp <- merge(heart4, liver4, by="SUBJID") 
length(unique(plotting_df_tmp$SUBJID)) # 148 individuals

# Split wide df into separate dfs so they can be recombined into a long format df
tmp1 <- plotting_df_tmp[,c("SUBJID", "gene.x", "organ.x", "SAMPID.x", "value.x")]
tmp2 <- plotting_df_tmp[,c("SUBJID", "gene.y", "organ.y", "SAMPID.y", "value.y")]

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

# Delete some old objects to open space
rm(voom_matrix)
rm(voom_df)
rm(tmp1)
rm(tmp2)

# Bind liver and heart dfs vertically 
plotting_df <- rbind(tmp1_2, tmp2_2)

#plotting_df <- merge(tmp1, tmp2, by=c("SUBJID"))
nrow(plotting_df)==nrow(tmp1_2)+nrow(tmp2_2) # TRUE
any(duplicated(plotting_df)) # FALSE

# Are the nrow as expected?
dim(plotting_df) # 15688 6
148 * 53 * 2 # 148 individuals x 53 genes x 2 samples per individual = 15688 expression values

# Transform the counts to log2(CPM)
plotting_df$log2_cpm <- cpm(plotting_df$value)

# Write plotting_df df to file
#write.table(plotting_df, "/scratch/mjpete11/linear_models/data/batch_voom_plotting_df.csv", sep=",")

# Read plotting_df df back in
#plotting_df <- read.csv("/scratch/mjpete11/linear_models/data/batch_voom_plotting_df.csv", sep=",")

# How many unique samples are there?
length(unique(plotting_df$SAMPID)) # 296

# Remove samples >6 standard deviations away
median(plotting_df$log2_cpm) # 78.61
mean(plotting_df$log2_cpm) # 63.74
sd(plotting_df$log2_cpm) # 77.09
sd(plotting_df$log2_cpm) * 6 # +/- 462.54
above <- 63.74 + 462.54 # 526.28
below <- 63.74 - 462.54 # -398.8
above;below

# Are there any samples outside of this range?
range(plotting_df$log2_cpm) # -146.90 234.91 

# Number of samples outside of this range
outlier_above <- plotting_df[plotting_df$log2_cpm > above,] # 0 samples 
outlier_below <- plotting_df[plotting_df$log2_cpm < below,] # 0 samples 
nrow(outlier_above);nrow(outlier_below)

######################### TEST PLOT ##########################################
dat <- plotting_df %>% filter(gene=="SLC25A1")
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
labs(x = "organ", y = "log2(CPM)", fill = "") +
ggtitle(paste0("Violin plot of ", GENE, " expression after blocking by batch via voom()")) 
ggsave(paste0("/scratch/mjpete11/linear_models/linear/voom_batch_violin_plots/", GENE, ".png"), device="png")


######################### END TEST ###########################################
# Function to plot violin plots
rm(plots)
rm(violin)
violin <- function(GENE){
		 dat <- plotting_df %>% filter(gene==GENE)
	     p <- ggplot(dat, aes(x = organ, y = log2_cpm, fill = organ)) +
				 	stat_compare_means(method = "wilcox.test",
									   aes(label = paste("adj.p_value =", after_stat(!!str2lang("p.adj"))*53)), 
									   label.x = 1.25, 
									   label.y = max(dat[["log2_cpm"]]) + 3,
									   paired = TRUE) +
					geom_violin(trim = FALSE) +
					stat_summary(fun.data = "mean_sdl", geom="crossbar", width=0.2, alpha=0.1) +
					scale_fill_manual(values = c("lightgreen", "purple")) +
					geom_jitter(size = 1, alpha = 0.9) +
					labs(x = "organ", y = "log2(CPM(prior.count=0.5))", fill = "") +
					ggtitle(paste0("Violin plot of ", GENE, " expression after blocking by batch via voom()")) 
			ggsave(paste0("/scratch/mjpete11/linear_models/batch_correction/voom_GTEx/voom_batch_violin_plots/", GENE, ".png"), device="png")
}
plots <- Map(violin, GENE=SLC)
plots
