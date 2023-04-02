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

# Delete some old objects to open space
rm(voom_mat)
rm(voom_matrix)
rm(tmp1)
rm(tmp2)

# Bind liver and heart dfs vertically 
organs <- rbind(tmp1_2, tmp2_2)

#organs <- merge(tmp1, tmp2, by=c("SUBJID"))
nrow(organs)==nrow(tmp1_2)+nrow(tmp2_2) # TRUE
any(duplicated(organs)) # FALSE

# Are the nrow as expected?
dim(organs) # 15688 5
148 * 53 * 2 # 148 individuals x 53 genes x 2 samples per individual = 15688 expression values

# Rename the 'value' column to 'log2_cpm'
colnames(organs)[5] <- "log2_cpm"

# Write organs df to file
write.table(organs, "/scratch/mjpete11/linear_models/data/batch_voom_organs.csv", sep=",")

# Read organs df back in
organs <- read.csv("/scratch/mjpete11/linear_models/data/batch_voom_organs.csv", sep=",")

# How many unique samples are there?
length(unique(organs$SAMPID)) # 296

# Write the striated samples to file
# Subset to the range of expected values
#striated <- subset(organs, organs$value > -7 & organs$value < -6.8)
striated <- organs[organs$log2_cpm < -6.8,]
unique(striated$gene) # UCP1, A52, A31, A2, A47, A48, A21, A41

# Subset the non-striated samples
'%ni%' <- Negate('%in%')
not_striated <- organs[organs$gene %ni% striated$gene,]

# Write distribution of gene values to file
write.table(subset(striated, gene=="UCP1"), "/scratch/mjpete11/linear_models/data/sparsity_matrix_files/UCP1_batch_voom.csv", sep=",", row.names=FALSE)
write.table(subset(striated, gene=="SLC25A52"), "/scratch/mjpete11/linear_models/data/sparsity_matrix_files/A52_batch_voom.csv", sep=",", row.names=FALSE)
write.table(subset(striated, gene=="SLC25A31"), "/scratch/mjpete11/linear_models/data/sparsity_matrix_files/A31_batch_voom.csv", sep=",", row.names=FALSE)
write.table(subset(striated, gene=="SLC25A2"), "/scratch/mjpete11/linear_models/data/sparsity_matrix_files/A2_batch_voom.csv", sep=",", row.names=FALSE)
write.table(subset(striated, gene=="SLC25A47"), "/scratch/mjpete11/linear_models/data/sparsity_matrix_files/A47_batch_voom.csv", sep=",", row.names=FALSE)
write.table(subset(striated, gene=="SLC25A48"), "/scratch/mjpete11/linear_models/data/sparsity_matrix_files/A48_batch_voom.csv", sep=",", row.names=FALSE)
write.table(subset(striated, gene=="SLC25A21"), "/scratch/mjpete11/linear_models/data/sparsity_matrix_files/A21_batch_voom.csv", sep=",", row.names=FALSE)
write.table(subset(striated, gene=="SLC25A41"), "/scratch/mjpete11/linear_models/data/sparsity_matrix_files/A41_batch_voom.csv", sep=",", row.names=FALSE)

# These are the genes from the non-batch corrected heart and liver samples
# that did not have an excess of zero values (i.e. the "non-striated" samples)
write.table(subset(not_striated, gene=="SLC25A1"), "/scratch/mjpete11/linear_models/data/sparsity_matrix_files/A1_batch.csv", sep=",", row.names=FALSE)
write.table(subset(not_striated, gene=="SLC25A3"), "/scratch/mjpete11/linear_models/data/sparsity_matrix_files/A3_batch.csv", sep=",", row.names=FALSE)
write.table(subset(not_striated, gene=="SLC25A4"), "/scratch/mjpete11/linear_models/data/sparsity_matrix_files/A4_batch.csv", sep=",", row.names=FALSE)
write.table(subset(not_striated, gene=="SLC25A5"), "/scratch/mjpete11/linear_models/data/sparsity_matrix_files/A5_batch.csv", sep=",", row.names=FALSE)
write.table(subset(not_striated, gene=="SLC25A6"), "/scratch/mjpete11/linear_models/data/sparsity_matrix_files/A6_batch.csv", sep=",", row.names=FALSE)
write.table(subset(not_striated, gene=="UCP2"), "/scratch/mjpete11/linear_models/data/sparsity_matrix_files/UCP2_batch.csv", sep=",", row.names=FALSE)
write.table(subset(not_striated, gene=="UCP3"), "/scratch/mjpete11/linear_models/data/sparsity_matrix_files/UCP3_batch.csv", sep=",", row.names=FALSE)
write.table(subset(not_striated, gene=="SLC25A10"), "/scratch/mjpete11/linear_models/data/sparsity_matrix_files/A10_batch.csv", sep=",", row.names=FALSE)
write.table(subset(not_striated, gene=="SLC25A11"), "/scratch/mjpete11/linear_models/data/sparsity_matrix_files/A11_batch.csv", sep=",", row.names=FALSE)
write.table(subset(not_striated, gene=="SLC25A12"), "/scratch/mjpete11/linear_models/data/sparsity_matrix_files/A12_batch.csv", sep=",", row.names=FALSE)
write.table(subset(not_striated, gene=="SLC25A13"), "/scratch/mjpete11/linear_models/data/sparsity_matrix_files/A13_batch.csv", sep=",", row.names=FALSE)
write.table(subset(not_striated, gene=="SLC25A14"), "/scratch/mjpete11/linear_models/data/sparsity_matrix_files/A14_batch.csv", sep=",", row.names=FALSE)
write.table(subset(not_striated, gene=="SLC25A15"), "/scratch/mjpete11/linear_models/data/sparsity_matrix_files/A15_batch.csv", sep=",", row.names=FALSE)
write.table(subset(not_striated, gene=="SLC25A16"), "/scratch/mjpete11/linear_models/data/sparsity_matrix_files/A16_batch.csv", sep=",", row.names=FALSE)
write.table(subset(not_striated, gene=="SLC25A17"), "/scratch/mjpete11/linear_models/data/sparsity_matrix_files/A17_batch.csv", sep=",", row.names=FALSE)
write.table(subset(not_striated, gene=="SLC25A18"), "/scratch/mjpete11/linear_models/data/sparsity_matrix_files/A18_batch.csv", sep=",", row.names=FALSE)
write.table(subset(not_striated, gene=="SLC25A19"), "/scratch/mjpete11/linear_models/data/sparsity_matrix_files/A19_batch.csv", sep=",", row.names=FALSE)
write.table(subset(not_striated, gene=="SLC25A20"), "/scratch/mjpete11/linear_models/data/sparsity_matrix_files/A20_batch.csv", sep=",", row.names=FALSE)
write.table(subset(not_striated, gene=="SLC25A22"), "/scratch/mjpete11/linear_models/data/sparsity_matrix_files/A22_batch.csv", sep=",", row.names=FALSE)
write.table(subset(not_striated, gene=="SLC25A23"), "/scratch/mjpete11/linear_models/data/sparsity_matrix_files/A23_batch.csv", sep=",", row.names=FALSE)
write.table(subset(not_striated, gene=="SLC25A24"), "/scratch/mjpete11/linear_models/data/sparsity_matrix_files/A24_batch.csv", sep=",", row.names=FALSE)
write.table(subset(not_striated, gene=="SLC25A25"), "/scratch/mjpete11/linear_models/data/sparsity_matrix_files/A25_batch.csv", sep=",", row.names=FALSE)
write.table(subset(not_striated, gene=="SLC25A26"), "/scratch/mjpete11/linear_models/data/sparsity_matrix_files/A26_batch.csv", sep=",", row.names=FALSE)
write.table(subset(not_striated, gene=="SLC25A27"), "/scratch/mjpete11/linear_models/data/sparsity_matrix_files/A27_batch.csv", sep=",", row.names=FALSE)
write.table(subset(not_striated, gene=="SLC25A28"), "/scratch/mjpete11/linear_models/data/sparsity_matrix_files/A28_batch.csv", sep=",", row.names=FALSE)
write.table(subset(not_striated, gene=="SLC25A29"), "/scratch/mjpete11/linear_models/data/sparsity_matrix_files/A29_batch.csv", sep=",", row.names=FALSE)
write.table(subset(not_striated, gene=="SLC25A30"), "/scratch/mjpete11/linear_models/data/sparsity_matrix_files/A30_batch.csv", sep=",", row.names=FALSE)
write.table(subset(not_striated, gene=="SLC25A32"), "/scratch/mjpete11/linear_models/data/sparsity_matrix_files/A32_batch.csv", sep=",", row.names=FALSE)
write.table(subset(not_striated, gene=="SLC25A33"), "/scratch/mjpete11/linear_models/data/sparsity_matrix_files/A33_batch.csv", sep=",", row.names=FALSE)
write.table(subset(not_striated, gene=="SLC25A34"), "/scratch/mjpete11/linear_models/data/sparsity_matrix_files/A34_batch.csv", sep=",", row.names=FALSE)
write.table(subset(not_striated, gene=="SLC25A35"), "/scratch/mjpete11/linear_models/data/sparsity_matrix_files/A35_batch.csv", sep=",", row.names=FALSE)
write.table(subset(not_striated, gene=="SLC25A36"), "/scratch/mjpete11/linear_models/data/sparsity_matrix_files/A36_batch.csv", sep=",", row.names=FALSE)
write.table(subset(not_striated, gene=="SLC25A37"), "/scratch/mjpete11/linear_models/data/sparsity_matrix_files/A37_batch.csv", sep=",", row.names=FALSE)
write.table(subset(not_striated, gene=="SLC25A38"), "/scratch/mjpete11/linear_models/data/sparsity_matrix_files/A38_batch.csv", sep=",", row.names=FALSE)
write.table(subset(not_striated, gene=="SLC25A39"), "/scratch/mjpete11/linear_models/data/sparsity_matrix_files/A39_batch.csv", sep=",", row.names=FALSE)
write.table(subset(not_striated, gene=="SLC25A40"), "/scratch/mjpete11/linear_models/data/sparsity_matrix_files/A40_batch.csv", sep=",", row.names=FALSE)
write.table(subset(not_striated, gene=="SLC25A42"), "/scratch/mjpete11/linear_models/data/sparsity_matrix_files/A42_batch.csv", sep=",", row.names=FALSE)
write.table(subset(not_striated, gene=="SLC25A43"), "/scratch/mjpete11/linear_models/data/sparsity_matrix_files/A43_batch.csv", sep=",", row.names=FALSE)
write.table(subset(not_striated, gene=="SLC25A44"), "/scratch/mjpete11/linear_models/data/sparsity_matrix_files/A44_batch.csv", sep=",", row.names=FALSE)
write.table(subset(not_striated, gene=="SLC25A45"), "/scratch/mjpete11/linear_models/data/sparsity_matrix_files/A45_batch.csv", sep=",", row.names=FALSE)
write.table(subset(not_striated, gene=="SLC25A46"), "/scratch/mjpete11/linear_models/data/sparsity_matrix_files/A46_batch.csv", sep=",", row.names=FALSE)
write.table(subset(not_striated, gene=="MTCH1"), "/scratch/mjpete11/linear_models/data/sparsity_matrix_files/MTCH1_batch.csv", sep=",", row.names=FALSE)
write.table(subset(not_striated, gene=="MTCH2"), "/scratch/mjpete11/linear_models/data/sparsity_matrix_files/MTCH2_batch.csv", sep=",", row.names=FALSE)
write.table(subset(not_striated, gene=="SLC25A49"), "/scratch/mjpete11/linear_models/data/sparsity_matrix_files/A49_batch.csv", sep=",", row.names=FALSE)
write.table(subset(not_striated, gene=="SLC25A50"), "/scratch/mjpete11/linear_models/data/sparsity_matrix_files/A50_batch.csv", sep=",", row.names=FALSE)
write.table(subset(not_striated, gene=="SLC25A51"), "/scratch/mjpete11/linear_models/data/sparsity_matrix_files/A51_batch.csv", sep=",", row.names=FALSE)
write.table(subset(not_striated, gene=="SLC25A53"), "/scratch/mjpete11/linear_models/data/sparsity_matrix_files/A53_batch.csv", sep=",", row.names=FALSE)

# Remove samples >6 standard deviations away
median(organs$value) # 3.74
mean(organs$value) # 3.04
sd(organs$value) # 3.67 
sd(organs$value) * 6 # +/- 22.04
above <- 2.85 + 23.3 # 26.15
below <- 2.85 - 23.3 # -20.45

# Are there any samples outside of this range?
range(organs$value) # -7.00 to 11.19 

# Number of samples outside of this range
outlier_above <- organs[organs$value > 26.15,] # 0 samples 
outlier_below <- organs[organs$value < -20.45,] # 0 samples 

# Subset one gene for the test plot
test <- subset(organs, gene %in% c("SLC25A6"))

# Test violin plot 
dat0 <- organs %>% filter(gene=="SLC25A6")
dat <- unique(dat0) 
p_val <- wilcox.test(formula=value~organ, data=dat, paired=TRUE, exact=TRUE)$p.value
n_tests <- 53
corrected_pval <- p.adjust(p_val, method="bonferroni", n=n_tests)
p <- ggplot(dat, aes(x = organ, y = value, fill = organ)) +
geom_violin(trim = FALSE) +
stat_summary(fun.data = "mean_sdl", geom="crossbar", width=0.2, alpha=0.1) +
scale_fill_manual(values = c("lightgreen", "purple")) +
geom_jitter(size = 1, alpha = 0.9) +
#scale_y_continuous(limits = c(-35, 35), expand = c(0,0),
#			       breaks = seq(-35, 35, by = 1)) +
labs(x = "organ", y = "log2(CPM)", fill = "") +
annotate(geom = "text", x = 1.5, y = 13, label=paste0("adj. p value: ", corrected_pval)) +
ggtitle(paste0("Violin plot of SLC25A6 expression between heart and liver after blocking by batch via limma::voom()")) 
ggsave(paste0("/scratch/mjpete11/linear_models/linear/voom_qnorm_violin_plots4/SLC25A6.png"), device="png")

# Function to plot violin plots
rm(plots)
rm(violin)
violin <- function(GENE){
		 dat <- organs %>% filter(gene==GENE)
         p_val <- wilcox.test(formula=log2_cpm~organ, data=dat, paired=TRUE, exact=TRUE)$p.value
	     n_tests <- 53
	     corrected_pval <- p.adjust(p_val, method="bonferroni", n=n_tests)
	     p <- ggplot(dat, aes(x = organ, y = log2_cpm, fill = organ)) +
	   	  	  geom_violin(trim = FALSE) +
       		  stat_summary(fun.data = "mean_sdl", geom="crossbar", width=0.2, alpha=0.1) +
	   		  scale_fill_manual(values = c("lightgreen", "purple")) +
	   		  geom_jitter(size = 1, alpha = 0.9) +
#       		  scale_y_continuous(limits = c(-35, 35), expand = c(0,0),
#			   			       breaks = seq(-35, 35, by = 1)) +
	   		  labs(x = "organ", y = "log2(CPM)", fill = "") +
			  annotate(geom = "text", x = 1.5, y = 13, label=paste0("adj. p value: ", corrected_pval)) +
	   		  ggtitle(paste0("Violin plot of ", GENE, " expression after blocking by batch via voom()")) 
	   		  ggsave(paste0("/scratch/mjpete11/linear_models/linear/voom_batch_violin_plots1/", GENE, ".png"), device="png")
}
plots <- Map(violin, GENE=SLC)
plots
