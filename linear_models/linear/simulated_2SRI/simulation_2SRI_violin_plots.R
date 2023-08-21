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

# Read in counts
counts <- fread("/scratch/mjpete11/linear_models/data/combined_simulated_001_batch1_batch2.csv") # integer
counts[1:5,1:5]

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

# Subset just 148 samples from each batch to match the GTEx analysis
test <- counts[,1:148+2]
test1 <- counts[,1003:1150]
head(colnames(test))
head(colnames(test))
tail(colnames(test1))
tail(colnames(test1))
test$ensembl_ID_no_version <- vec
test1$ensembl_ID_no_version <- vec
counts_subset <- merge(test, test1, by=c("ensembl_ID_no_version"))
counts_subset[1:5,1:5]
dim(counts_subset) # 61386 297

# Subset the SLC25 genes
sub_df <- counts_subset[counts_subset$"ensembl_ID_no_version" %in% SLC25, ]

# Is every gene present? 
setdiff(SLC25, sub_df$'ensembl_ID_no_version') # Victory!

# Number of genes remaining
nrow(sub_df) # 53 

# Reshape from long to wide format
# Add gene ID column
sub_df$hugo_ID <- SLC
sub_df <- sub_df %>% select(hugo_ID, everything())
sub_df$ensembl_ID <- NULL
sub_df[1:5,1:5]

# Reshape from wide to long format so the dataframe is compartable with the 
# plotting function
# Reshape dataframe so it can be converted to a design matrix object
plotting_df <- sub_df %>% reshape2::melt(id.vars=c("hugo_ID","ensembl_ID_no_version"))
head(plotting_df)
dim(plotting_df) # 15582 4

# Rename column labeling the batch
colnames(plotting_df)[3] <- c("sample")

# Add a column with the batch
batch1 <- rep_len('batch1', len=nrow(plotting_df)/2)
batch2 <- rep_len('batch2', len=nrow(plotting_df)/2)
plotting_df$batch <- c(batch1, batch2)
head(plotting_df)
tail(plotting_df)

# Add a column with the log2(CPM) transformed values
plotting_df$log2_cpm <- cpm(plotting_df$value)
head(plotting_df)

# Remove samples >6 standard deviations away
median(plotting_df$log2_cpm) # 54.0
mean(plotting_df$log2_cpm) # 64.18
sd(plotting_df$log2_cpm) # 33.34
sd(plotting_df$log2_cpm) * 6 # +/- 200.02 
above <- 3.43 + 20.6 # 24.03
below <- 3.43 - 20.6 # -17.17
above;below

# Are there any samples outside of this range?
range(plotting_df$log2_cpm) # 0.00 to 176.07 

# Number of samples outside of this range
outlier_above <- plotting_df[plotting_df$log2_cpm > 24.03,] # 0 sample  
outlier_below <- plotting_df[plotting_df$log2_cpm < -17.77,] # 0 samples 
nrow(outlier_below);nrow(outlier_below)

# Make design matrix; want to end up with 148 paired heart and liver samples
# to match the GTEx design
# Subset just 148 samples from each batch to match the GTEx analysis
counts_subset[1:5,1:5]
counts_subset$ensembl_ID_no_version <- NULL  # Drop the ensembl_ID column before applying voom
# since you need a numeric matrix only
dim(counts_subset) # 61386 296
counts_subset[1:5,1:5]
heart_vec <- rep_len('heart', len=148)
liver_vec <- rep_len('liver', len=148)
length(heart_vec) # 148 
length(liver_vec) # 148
organ_vect <- c(heart_vec, liver_vec)
mat <- fastDummies::dummy_cols(organ_vect)
head(mat);tail(mat)
class(mat) # data.frame
mat <- data.matrix(mat)
class(mat) # matrix
dim(mat) # 592 3
nrow(mat);ncol(counts_subset) # 296 296
nrow(mat)==ncol(counts_subset) # TRUE

################################### TEST ###################################
# Create test counts and design matrix
#test_counts <- counts[1:10,1:10]
# Add some non-zero values
#test_counts[1:5,1] <- sample.int(1000, 5, 3)
#test_counts[1:5,9] <- sample.int(1000, 5, 3)
#test_counts[6:10,3] <- sample.int(1000, 5, 3)
#test_counts[6:10,6] <- sample.int(1000, 5, 3)
#test_counts

# Change the last 5 rows to indicate liver instead of heart samples
# otherwise voom will not run since two conditions are not present
#test_mat <- mat[1:10,]
#test_mat[6:10,3] <- 1 
#test_mat[6:10,2] <- 0
#test_mat

# Apply limma::voom() to counts matrix to normalize and write to file
#voom_obj <- voom(counts=test_counts, design=test_mat, normalize.method="quantile")  
#write.csv("/scratch/mjpete11/mitochondrial_paralogs/data/data/test_simulated_data_voom_Eobject.tsv", sep="\t") # float
#print("completed voom")
#head(voom_obj$E)
################################### TEST ###################################

# Apply limma::voom() to counts matrix to normalize and write to file
voom_obj <- voom(counts=counts_subset, design=mat, normalize.method="quantile")  
head(voom_obj$E)

# Store voom adjusted values
voom_E <- as.data.frame(voom_obj$E)
voom_E$ensembl_IDs  <- vec 
voom_E[1:5,1:5]
# Move the ensembl <- ID column to the front
voom_E <- voom_E %>%
		  dplyr::select(ensembl_IDs, everything())
voom_E[1:5,1:5]
dim(voom_E) # 61386 297
class(voom_E)
#write.csv(voom_E, "/scratch/mjpete11/linear_models/data/simulated_data_voom_Eobject.csv") # float
print("completed voom")

# Subset the gene_counts df using the ensembl_IDs; then use the corresponding
# index to subset the voom-transformed matrix (plotting_df)
gene_counts <- voom_E[voom_E$ensembl_IDs %in% SLC25,]
gene_counts[1:5,1:5]
dim(gene_counts) # 53 297

# Add gene ID column
gene_counts$hugo_ID <- SLC
gene_counts[1:5,1:5]

# Reshape from wide to long format so the dataframe is compartable with the 
# plotting function
# Reshape dataframe so it can be converted to a design matrix object
plotting_df <- gene_counts %>% reshape2::melt(id.vars=c("hugo_ID","ensembl_IDs"))
head(plotting_df)
dim(plotting_df) # 15688 5

# Add a column with the 'organ' status at random but with equal heart and liver
# sample sizes
heart_vec <- rep_len('heart', len=nrow(plotting_df)/2)
liver_vec <- rep_len('liver', len=nrow(plotting_df)/2)
length(heart_vec) # 7844 
length(liver_vec) # 7844
plotting_df$organ <- c(heart_vec, liver_vec)
head(plotting_df)
tail(plotting_df)

# Simulate a range of RIN and ischemic time values and add to dataframe
#Read in voom quantile normalized counts df with GTEx RIN and ischemic time values
organs <- fread("/scratch/mjpete11/linear_models/data/voom_qnorm_counts1.csv", sep=",") # float
range(organs$SMRIN) # 5.6, 9.6
range(organs$SMTSISCH) # 69, 1426
rin <- round(runif(nrow(plotting_df), 5.6, 9.6), digits=1) 
ischemic_time <- round(runif(nrow(plotting_df), 69, 1426))
plotting_df$SMRIN <- rin
plotting_df$SMTSISCH <- ischemic_time
head(plotting_df)
tail(plotting_df)

# Add column with the batch label
# Add a column with the batch
batch1 <- rep_len('batch1', len=nrow(plotting_df)/2)
batch2 <- rep_len('batch2', len=nrow(plotting_df)/2)
plotting_df$batch <- c(batch1, batch2)
head(plotting_df)
tail(plotting_df)

# Add a column with the log2(CPM) transformed values
plotting_df$log2_cpm <- cpm(plotting_df$value)
head(plotting_df)

# Set y lim
range(plotting_df$value) # c(4.10, 15.47) at 5 counts filtering threshold

# before making the violin plots
colnames(plotting_df) # "hugo_ID", "ensembl_IDs", "variable", "value", "organ", "SMRIN", "SMTSISCH"

##### TEST; Make one dataframe #####
#datt <- plotting_df %>% filter(hugo_ID == "SLC25A1")
#datt2 <- na.omit(datt)
# Regress out effect of RIN and ischemic time on SLC
#model1 <- lm(value ~ SMRIN + SMTSISCH, data=datt2)
# Add fitted values to dataframe
#datt2$resids <- residuals(model1)
# Regress organ on SLC using previous fitted regression values as offset
#model2 <- lm(value ~ batch + resids, data=datt2)
# Add fitted values from second linear model to dattaframs
#datt2$fitted_values <- fitted.values(model2)
#head(datt2)
##### TEST #####

# Make a list of dataframes with the 2SRI corrected values
contained <- list()
for (paralog in SLC){
	datt <- plotting_df %>% filter(hugo_ID == paralog)
	datt2 <- na.omit(datt)
	datt3 <- distinct(na.omit(datt2))
	# Regress out effect of RIN and ischemic time on SLC
	model1 <- lm(value ~ SMRIN + SMTSISCH, data=datt3)
	# Add fitted values to dataframe
	datt3$resids <- residuals(model1)
	# Regress organ on SLC using previous fitted regression values as offset
	model2 <- lm(value ~ batch + resids, data=datt3)
	# Add fitted values from second linear model to data frames
	datt3$fitted_values <- fitted.values(model2)
	contained[[paralog]] <- datt3
}
organs1 <- ldply(contained, data.frame)  
head(organs1)
tail(organs1)

# Remove samples >6 standard deviations away
median(organs1$fitted_values) # 13.95
mean(organs1$fitted_values) # 13.88
sd(organs1$fitted_values) # 1.51
sd(organs1$fitted_values) * 6 # +/- 9.07 
above <- 3.96 + 24.0 # 27.96
below <- 3.96 - 24.0 # -20.04
range(organs1$fitted_values) # 4.09 14.5

# Violin and jitter plot function
rm(violin_plots)

violin_plots <- function(GENE){
	dat <- organs1 %>% filter(hugo_ID == GENE)
    dat2 <- na.omit(dat)
    dat3 <- distinct(na.omit(dat2))
    # Regress out effect of RIN and ischemic time on SLC
    model1 <- lm(value ~ SMRIN + SMTSISCH, data=dat3)
	# Add fitted values to dataframe
	dat3$resids <- residuals(model1)
	# Regress organ on SLC using previous fitted regression values as offset
	model2 <- lm(value ~ batch + resids, data=dat3)
	# Add fitted values from second linear model to dataframs
	dat3$fitted_values <- fitted.values(model2)
	p_val <- wilcox.test(formula=fitted_values~batch,data=dat3)$p.value
	corrected_pval <- p.adjust(p_val, method="bonferroni", n=53)
    # Violin plot
	p <- ggplot(dat3, aes(x = batch, y = fitted_values, fill=batch)) +
	geom_violin(width=1, position=position_dodge(width=0.5)) +
	scale_fill_manual(values=c("lightgreen", "purple")) +
	geom_jitter(width=0.3) +
	stat_summary(fun.data="mean_sdl", geom="crossbar", width=0.1, alpha=0.1) +
	#annotate(geom = "text", x = 1.5, y = 19, label=paste0("p value: ",corrected_pval)) +
	annotate(geom = "text", x = 1.5, y = 34, label=paste0("p value: ",corrected_pval)) +
	#ylim(c(0,30)) +
	ylab("log2(CPM)") +
	xlab("organ") +
#	scale_y_continuous(limits = c(-35, 35), expand = c(0,0), breaks = seq(-35, 35, by = 1)) +
#	scale_y_continuous(breaks=scales::pretty_breaks(n=20)) + # Add more y-axis tick marks
	ggtitle(paste0("Violin plot of simulated ", GENE, "\n expression between heart and liver after 2SRI")) +
	theme(plot.title = element_text(hjust = 0.5, size=14)) 
	ggsave(paste0("/scratch/mjpete11/linear_models/linear/simulated_2SRI/plots_error_001/", GENE, ".png"), p, device="png")
}
plots <- Map(violin_plots, GENE=SLC)
violin_plots(GENE="SLC25A1", DATA=organs)
plots[[1]]
plots
