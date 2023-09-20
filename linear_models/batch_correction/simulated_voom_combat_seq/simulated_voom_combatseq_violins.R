# Purpose: Linear models to assess if SLC paralogs are dependent on energetic state?

# Libraries
library(plyr)
library(dplyr)
library(data.table)
library(tidyverse)
library(reshape)
library(stringr)
library(ggplot2)
library(ggpubr)
library(rstatix)
library(stats)
library(scales)
library(limma)
library(fastDummies)
library(edgeR)
library(sva)

# Read in counts
counts <- fread("/scratch/mjpete11/linear_models/data/combined_simulated_001_batch1_batch2.csv") # integer
counts[1:5,1:5]
dim(counts) # 61386 2001

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

SLC25_with_version <- c("ENSG00000100075.10", "ENSG00000120329.7", "ENSG00000075415.15",
		   "ENSG00000151729.11", "ENSG00000005022.6", "ENSG00000169100.14",
		   "ENSG00000109424.4", "ENSG00000175567.11", "ENSG00000175564.13",
		   "ENSG00000183048.12", "ENSG00000108528.14", "ENSG00000115840.14",
		   "ENSG00000004864.14", "ENSG00000102078.16", "ENSG00000102743.15",
		   "ENSG00000122912.15", "ENSG00000100372.15", "ENSG00000182902.14",
		   "ENSG00000125454.12", "ENSG00000178537.10", "ENSG00000183032.12",
		   "ENSG00000177542.11", "ENSG00000125648.15", "ENSG00000085491.17",
		   "ENSG00000148339.13", "ENSG00000144741.19", "ENSG00000153291.16",
		   "ENSG00000155287.11", "ENSG00000197119.13", "ENSG00000174032.17",
		   "ENSG00000151475.6", "ENSG00000164933.12", "ENSG00000171612.7",
		   "ENSG00000162461.8", "ENSG00000125434.11", "ENSG00000114120.14",
		   "ENSG00000147454.14", "ENSG00000144659.13", "ENSG00000013306.16",
		   "ENSG00000075303.13", "ENSG00000181240.14", "ENSG00000181035.14",
		   "ENSG00000077713.19", "ENSG00000160785.14", "ENSG00000162241.13",
		   "ENSG00000164209.17", "ENSG00000140107.12", "ENSG00000145832.15",
		   "ENSG00000137409.20", "ENSG00000109919.10", "ENSG00000122696.14",
		   "ENSG00000141437.10", "ENSG00000269743.3")

SLC25 <- c("ENSG00000100075", "ENSG00000120329", "ENSG00000075415",
		   "ENSG00000151729", "ENSG00000005022", "ENSG00000169100",
		   "ENSG00000109424", "ENSG00000175567", "ENSG00000175564",
		   "ENSG00000183048", "ENSG00000108528", "ENSG00000115840",
		   "ENSG00000004864", "ENSG00000102078", "ENSG00000102743",
		   "ENSG00000122912", "ENSG00000100372", "ENSG00000182902",
		   "ENSG00000125454", "ENSG00000178537", "ENSG00000183032",
		   "ENSG00000177542", "ENSG00000125648", "ENSG00000085491",
		   "ENSG00000148339", "ENSG00000144741", "ENSG00000153291",
		   "ENSG00000155287", "ENSG00000197119", "ENSG00000174032",
		   "ENSG00000151475", "ENSG00000164933", "ENSG00000171612",
		   "ENSG00000162461", "ENSG00000125434", "ENSG00000114120",
		   "ENSG00000147454", "ENSG00000144659", "ENSG00000013306",
		   "ENSG00000075303", "ENSG00000181240", "ENSG00000181035",
		   "ENSG00000077713", "ENSG00000160785", "ENSG00000162241",
		   "ENSG00000164209", "ENSG00000140107", "ENSG00000145832",
		   "ENSG00000137409", "ENSG00000109919", "ENSG00000122696",
		   "ENSG00000141437", "ENSG00000269743")

# Convert to df instead of df + datatable
class(counts)
counts <- as.data.frame(counts) 

# Drop the version number (decimals) after the ensembl_ID
counts$ensembl_ID <- as.character(counts$ensembl_ID) # Convert from factor to character 
class(counts$ensembl_ID) # character
vec <- gsub("\\.\\d+", "", counts$ensembl_ID) # write this to file because it takes forever to run
#write.csv(vec, "/scratch/mjpete11/linear_models/data/ensembl_IDs_no_version.csv")
#vec <- read.csv("/scratch/mjpete11/linear_models/data/ensembl_IDs_no_version.csv")
counts$ensembl_ID_no_version <- vec
counts[1:5,1:5]

# Move last column first
counts <- counts %>% select(ensembl_ID_no_version, everything())
counts[1:5,1:5]
dim(counts) # 61386 2002

# Subset just 148 samples from each batch to match the GTEx analysis
test <- counts[,1:148+2]
test1 <- counts[,1003:1150]
head(colnames(test))
head(colnames(test1))
tail(colnames(test))
tail(colnames(test1))
test$ensembl_ID_no_version <- vec
test1$ensembl_ID_no_version <- vec
counts_subset <- merge(test, test1, by=c("ensembl_ID_no_version"))
counts_subset[1:5,1:5]
dim(counts_subset) # 61386 297
head(colnames(counts_subset))
tail(colnames(counts_subset))

# Make design matrix; want to end up with 148 paired heart and liver samples
# to match the GTEx design
# Subset just 148 samples from each batch to match the GTEx analysis
counts_subset[1:5,1:5]
counts_subset$ensembl_ID_no_version <- NULL  # Drop the ensembl_ID column before applying voom
# since you need a numeric matrix only
dim(counts_subset) # 61386 296
counts_subset[1:5,1:5]

# Make design matrix
batch1_vec <- rep_len('batch1', len=148)
batch2_vec <- rep_len('batch2', len=148)
length(batch1_vec) # 148 
length(batch2_vec) # 148
batch_vect <- c(batch1_vec, batch2_vec)
design_matrix <- as.data.frame(model.matrix(~batch_vect-1))
# Rename columns so they are more clear
colnames(design_matrix) <- c("batch1", "batch2")
nrow(design_matrix)==ncol(counts_subset) # TRUE
dim(counts_subset) # 61386 296

# Apply limma::voom() to counts matrix to normalize and write to file
voom_obj <- voom(counts=counts_subset, design=design_matrix, normalize.method="quantile")  
head(voom_obj$E)

# Store voom adjusted values
voom_E <- as.data.frame(voom_obj$E)
voom_E$ensembl_IDs  <- vec 
voom_E[1:5,1:5]
# Move the ensembl_ID column to the front
voom_E <- voom_E %>%
		  dplyr::select(ensembl_IDs, everything())
voom_E[1:5,1:5]
dim(voom_E) # 61386 297
class(voom_E)
#write.csv(voom_E, "/scratch/mjpete11/linear_models/data/simulated_error_5_voom_E_object.csv") # float
print("completed voom")

# Make design vector for the filtered samples dataframe
batch1_vec <- rep_len('batch1', len=ncol(voom_E)/2)
batch2_vec <- rep_len('batch2', len=ncol(voom_E)/2)
length(batch1_vec) # 148 
length(batch2_vec) # 148
batch_vect <- c(batch1_vec, batch2_vec)
length(batch_vect)==ncol(voom_E)-1 # TRUE
sapply(batch_vect, class)
batch_vect <- as.factor(batch_vect)
sapply(batch_vect, class)
head(batch_vect)
batch_dummy_vec <- as.numeric(batch_vect)
head(batch_dummy_vec);tail(batch_dummy_vec)

# Apply ComBat_seq to the voom adjusted data, using parametric empirical Bayesian adjustment
combat_edata <- ComBat_seq(counts=voom_obj, batch=batch_dummy_vec)
class(combat_edata) # matrix

# Append the ensembl_IDs on the dataframe
combat_edata <- as.data.frame(combat_edata)
combat_edata$ensembl_IDs  <- vec 
combat_edata[1:5,1:5]
# Move the ensembl_ID column to the front
combat_edata <- combat_edata %>%
		  dplyr::select(ensembl_IDs, everything())
combat_edata[1:5,1:5]
dim(combat_edata) # 61386 297
class(combat_edata) # data.frame

# Subset the gene_counts df using the ensembl_IDs; then use the corresponding
# index to subset the voom-transformed matrix (plotting_df)
gene_counts <- combat_edata[combat_edata$ensembl_IDs %in% SLC25,]
gene_counts[1:5,1:5]
dim(gene_counts) # 53 297

# Add gene ID column
gene_counts$hugo_ID <- SLC
gene_counts <- gene_counts %>%
		  dplyr::select(hugo_ID, everything())
gene_counts[1:5,1:5]

# Reshape from wide to long format so the dataframe is compartable with the 
# plotting function
# Reshape dataframe so it can be converted to a design matrix object
plotting_df <- gene_counts %>% reshape2::melt(id.vars=c("hugo_ID","ensembl_IDs"))
head(plotting_df)
dim(plotting_df) # 15688 4

# Rename column labeling the batch
colnames(plotting_df)[3] <- c("sample")

# Add a column with the batch to the plotting dataframe
batch1_vect <- rep_len('batch1', len=nrow(plotting_df)/2)
batch2_vect <- rep_len('batch2', len=nrow(plotting_df)/2)
length(batch1_vect) # 7844 
length(batch2_vect) # 7844
batch_vec <- c(batch1_vect, batch2_vect)
plotting_df$batch <- batch_vec
head(plotting_df);tail(plotting_df)

# Add a column with the log2(CPM) transformed value
#### skipping this step because the voom() function inclused the cpm() function
#plotting_df$log2_cpm <- cpm(plotting_df$value, log=TRUE, prior.count=0.5)
#head(plotting_df)

# Rename the 'value' column to log2_cpm
colnames(plotting_df)[4] <- "log2_cpm"

# Remove samples >6 standard deviations away
median(plotting_df$log2_cpm) # 14 
mean(plotting_df$log2_cpm) # 13.91
sd(plotting_df$log2_cpm) # 1.50
sd(plotting_df$log2_cpm) * 6 # +/- 9.01   
above <- mean(plotting_df$log2_cpm) + sd(plotting_df$log2_cpm) * 6
below <- mean(plotting_df$log2_cpm) - sd(plotting_df$log2_cpm) * 6
above;below # 22.92, 4.90

# Are there any samples outside of this range?
range(plotting_df$log2_cpm) # 5 16

# Number of samples outside of this range
outlier_above <- plotting_df[plotting_df$log2_cpm > above,] # 0 sample  
outlier_below <- plotting_df[plotting_df$log2_cpm < below,] # 0 samples 
nrow(outlier_below);nrow(outlier_below)

# Convert 'batch' variable to factor and 'log2_cpm' to numeric
# This is probably driving the errors with the p-value calulcation in
# stat_compare_means()
lapply(plotting_df$batch, class) # character
lapply(plotting_df$log2_cpm, class) # numeric

#plotting_df$batch <- lapply(plotting_df$batch, as.factor) # character
#lapply(plotting_df$batch, class) # factor 
# converting to factor causes: Error: Discrete value supplied to continuous scale

# Function to plot violin plots
rm(plots)
rm(violin)

violin <- function(GENE){
		dat <- plotting_df %>% filter(hugo_ID==GENE)
		p <- ggplot(dat, aes(x = batch, y = log2_cpm, fill = batch)) +
				stat_compare_means(method = "wilcox.test", 
								   aes(label = paste("adj.p_value =", after_stat(!!str2lang("p.adj"))*53)), 
								   label.x = 1.25, 
								   label.y = max(dat[["log2_cpm"]]) + 10,
								   label.y.npc = "top",
								   inherit.aes = TRUE,
								   paired = TRUE) +
				geom_violin(trim = FALSE) +
				stat_summary(fun.data = "mean_sdl", geom="crossbar", width=0.2, alpha=0.1) +
				scale_fill_manual(values = c("lightgreen", "purple")) +
				geom_jitter(size = 1, alpha = 0.9) +
				labs(x = "batch", y = "log2(CPM(prior.count=0.05))", fill = "") +
				scale_x_discrete(labels = c("batch 1: error_rate = 0.005", "batch 2: error_rate = 0.001")) +
				ggtitle(paste0("Violin plot of simulated ", GENE, "\n expression after adjustment via voom and combat_seq")) +
				theme(plot.title=element_text(hjust=0.5))
		ggsave(paste0("/scratch/mjpete11/linear_models/batch_correction/simulated_voom_combat_seq/plots_error_001/", GENE, ".png"), device="png")
}
plots <- Map(violin, GENE=SLC)

############################# TEST PLOT ######################################
dat <- plotting_df %>% filter(hugo_ID=="SLC25A2")
p <- ggplot(dat, aes(x = batch, y = log2_cpm, fill = batch)) +
		stat_compare_means(method = "wilcox.test", 
						   aes(label = paste("adj.p_value =", after_stat(!!str2lang("p.adj"))*53)), 
						   label.x = 1.25, 
						   label.y = max(dat[["log2_cpm"]]) + 0.5,
						   paired = TRUE) +
		geom_violin(trim = FALSE) +
		stat_summary(fun.data = "mean_sdl", geom="crossbar", width=0.2, alpha=0.1) +
		scale_fill_manual(values = c("lightgreen", "purple")) +
		geom_jitter(size = 1, alpha = 0.9) +
		labs(x = "batch", y = "log2(CPM(prior.count=0.5))", fill = "") +
		scale_x_discrete(labels = c("batch 1: error_rate = 0.005", "batch 2: error_rate = 0.001")) +
		ggtitle(paste0("Violin plot of simulated ", "SLC25A1", "\n expression after adjustment via voom and combat_seq")) +
		theme(plot.title=element_text(hjust=0.5))
ggsave(paste0("/scratch/mjpete11/linear_models/batch_correction/simulated_voom_combat_seq/plots_error_001/", "SLC25A2", ".png"), device="png")

rm(list=ls("p", "dat"))
############################# TEST PLOT ######################################
