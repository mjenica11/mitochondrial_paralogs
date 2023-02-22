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

# Read in voom + qnorm adjusted filtered GTEx count matrix
voom_matrix <- fread("/scratch/mjpete11/linear_models/data/batch_voom_qnorm_matrix.csv", sep=",") 

# Read in sample attributes files
file2 <- read.csv("/scratch/mjpete11/linear_models/data/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt", sep="\t")

# Add gene name column
gene_names <- fread("/scratch/mjpete11/linear_models/data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct", sep="\t",
					select=c("Description"))
gene_names <- as.vector(gene_names[,1])
voom_matrix$"Description" <- gene_names 

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

# Subset the SLC25 genes
voom_matrix <- as.data.frame(voom_matrix)
sub_df <- voom_matrix[voom_matrix$"Description" %in% SLC, ]

# Number of genes remaining
dim(voom_matrix) # 56200 15985
dim(sub_df) # 53 15985 

# Keep only unique colnames
voom_mat <- subset(sub_df, select=unique(colnames(voom_matrix)))
dim(voom_mat) # 53 297

# Move the gene names column to the front
voom_mat <- voom_mat %>% select("Description", everything())
voom_mat[1:5,1:5]

# Subset heart left ventricle and liver samples into two separate dfs
heart <- file2[file2$SMTSD %in% "Heart - Left Ventricle", ]
liver <- file2[file2$SMTSD %in% "Liver", ]

# Subset the voom + qnorm adjusted matrix by the sample IDs that match the IDs in file3 df
heart2 <- voom_mat %>% select(contains(heart$SAMPID))
liver2 <- voom_mat %>% select(contains(liver$SAMPID))

# How many heart/liver samples are there?
ncol(heart2) # 148
ncol(liver2) # 148

# Append the gene ensemble ID and common name
heart3 <- cbind('gene'=voom_mat$'Description', heart2)
liver3 <- cbind('gene'=voom_mat$'Description', liver2)

# Add a columns with the organ in each df
heart3$organ <- 'heart' 
liver3$organ <- 'liver' 

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
organs_tmp <- merge(heart4, liver4, by="SUBJID") 
length(unique(organs_tmp$SUBJID)) # 148 individuals
tmp1 <- organs_tmp[,c("SUBJID", "gene.x", "organ.x", "SAMPID.x", "value.x")]
tmp2 <- organs_tmp[,c("SUBJID", "gene.y", "organ.y", "SAMPID.y", "value.y")]
colnames(tmp1) <- c("SUBJID", "gene", "organ", "SUBJID", "value")
colnames(tmp2) <- c("SUBJID", "gene", "organ", "SUBJID", "value")
organs0 <- rbind(tmp1, tmp2)
colnames(organs0)[4] <- "SAMPID"
organs <- organs0[!duplicated(organs0),]

# Are the nrow as expected?
dim(organs) # 15688 5
148 * 53 * 2 # 148 individuals x 53 genes x 2 samples per individual = 15688 expression values

# Write organs df to file
write.table(organs, "/scratch/mjpete11/linear_models/data/batch_voom_organs.csv", sep=",")

# Read organs df back in
organs <- read.csv("/scratch/mjpete11/linear_models/data/batch_voom_organs.csv", sep=",")

# Write the striated samples to file
# Subset to the range of expected values
striated <- subset(organs, organs$value > -7 & organs$value < -6.8)
nrow(striated) # 967 genes
# Write to file
write.table(striated, "/scratch/mjpete11/linear_models/data/striated_voom_batch.csv", sep=",")

# Remove samples >6 standard deviations away
median(organs$value) # 3.68
sd(organs$value) # 3.88 
sd(organs$value) * 6 # +/- 23.3

# Are there any samples outside of this range?
range(organs$value) # -7.00 to 11.19 

# Number of samples outside of this range
outlier_above <- organs[organs$value > 23.3,] # 0 samples 
outlier_below <- organs[organs$value < -23.3,] # 0 samples 

# Subset one gene for the test plot
test <- subset(organs, gene %in% c("SLC25A4"))

# Test violin plot with boxplots of standard deviations
ggplot(test, aes(x = organ, y = value, fill = organ)) +
	   geom_violin(trim = FALSE) +
       stat_summary(fun.data = "mean_sdl", geom="crossbar", width=0.2) +
	   scale_fill_manual(values = c("lightgreen", "purple")) +
	   geom_jitter(size = 3, position = position_dodge(0.7), alpha = 0.3) +
       scale_y_continuous(limits = c(3, 12), expand = c(0,0),
						  breaks = round(seq(min(organs$value), max(organs$value), by = 0.5))) +
	   labs(x = "organ", y = "log2(CPM)", fill = "")
	   ggtitle(paste0("Violin plot of SLC expression between heart and liver after blocking by batch via limma::voom()")) 
	   ggsave(paste0("/scratch/mjpete11/linear_models/linear/voom_qnorm_violin_plots4/test.png"), device="png")

# Function to plot violin plots
rm(violin)
violin <- function(GENE){
		 dat <- organs %>% filter(gene==GENE)
         p_val <- wilcox.test(formula=value~organ, data=dat, paired=TRUE, exact=TRUE)$p.value
	     n_tests <- 53
	     corrected_pval <- p.adjust(p_val, method="bonferroni", n=n_tests)
	     p <- ggplot(dat, aes(x = organ, y = value, fill = organ)) +
	   	  	  geom_violin(trim = FALSE) +
       		  stat_summary(fun.data = "mean_sdl", geom="crossbar", width=0.2, alpha=0.1) +
	   		  scale_fill_manual(values = c("lightgreen", "purple")) +
	   		  geom_jitter(size = 1, alpha = 0.9) +
#       		  scale_y_continuous(limits = c(-35, 35), expand = c(0,0),
#			   			       breaks = seq(-35, 35, by = 1)) +
	   		  labs(x = "organ", y = "log2(CPM)", fill = "") +
			  annotate(geom = "text", x = 1.5, y = 13, label=paste0("adj. p value: ", corrected_pval)) +
	   		  ggtitle(paste0("Violin plot of ", GENE, " expression after blocking by batch via voom()")) 
	   		  ggsave(paste0("/scratch/mjpete11/linear_models/linear/voom_qnorm_violin_plots2/", GENE, ".png"), device="png")
}
plots <- Map(violin, GENE=SLC)
plots
