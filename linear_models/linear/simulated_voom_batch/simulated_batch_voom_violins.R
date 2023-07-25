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
voom_E <- fread("/scratch/mjpete11/linear_models/data/simulated_data_voom_Eobject.csv", sep=",") 
voom_E[1:10,1:10]

# Ensembl gene IDs
SLC25_IDs <- c("ENSG00000100075", "ENSG00000120329", "ENSG00000075415", "ENSG00000151729", "ENSG00000005022", "ENSG00000169100", "ENSG00000109424","ENSG00000175567", "ENSG00000175564", "ENSG00000183048", "ENSG00000108528", "ENSG00000115840", "ENSG00000004864", "ENSG00000102078", "ENSG00000102743", "ENSG00000122912", "ENSG00000100372", "ENSG00000182902", "ENSG00000125454", "ENSG00000178537", "ENSG00000183032","ENSG00000177542", "ENSG00000125648", "ENSG00000085491", "ENSG00000148339", "ENSG00000144741", "ENSG00000153291", "ENSG00000155287", "ENSG00000197119", "ENSG00000174032", "ENSG00000151475", "ENSG00000164933", "ENSG00000171612", "ENSG00000162461", "ENSG00000125434", "ENSG00000114120", "ENSG00000147454", "ENSG00000144659", "ENSG00000013306", "ENSG00000075303", "ENSG00000181240", "ENSG00000181035","ENSG00000077713", "ENSG00000160785", "ENSG00000162241", "ENSG00000164209", "ENSG00000140107", "ENSG00000145832", "ENSG00000137409","ENSG00000109919", "ENSG00000122696", "ENSG00000141437", "ENSG00000269743")

# Subset the gene_counts df using the ensembl_IDs; then use the corresponding
# index to subset the voom-transformed matrix (plotting_df)
gene_counts <- voom_E[voom_E$ensembl_IDs %in% SLC25_IDs,]
gene_counts[1:5,1:5]
dim(gene_counts) # 53 297

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

# Set y lim
range(plotting_df$value) # c(4.10, 15.47) at 5 counts filtering threshold

# before making the violin plots
colnames(plotting_df) # "hugo_ID", "ensembl_IDs", "variable", "value", "organ", "SMRIN", "SMTSISCH"
head(plotting_df)

# Remove samples >6 standard deviations away
median(plotting_df$value) # 13.95
mean(plotting_df$value) # 13.88
sd(plotting_df$value) # 1.51 
sd(plotting_df$value) * 6 # +/- 9.07
above <- 1.51 + 9.07 # 10.58 
below <- 1.51 - 9.07 # -7.56

# Are there any samples outside of this range?
range(plotting_df$value) # 4.09 15.47 

# Number of samples outside of this range
outlier_above <- plotting_df[plotting_df$value > 10.58,] # 0 samples 
outlier_below <- plotting_df[plotting_df$value < -7.56,] # 0 samples 

# Subset one gene for the test plot
test <- subset(plotting_df, hugo_ID %in% c("SLC25A6"))

# Test violin plot 
dat0 <- plotting_df %>% filter(hugo_ID=="SLC25A6")
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
ggtitle(paste0("Violin plot of simulated SLC25A6 expression between /n heart and liver after blocking by batch via limma::voom()")) 
ggsave(paste0("/scratch/mjpete11/linear_models/linear/simulated_voom_blocking_with_batch_plots/SLC25A6.png"), device="png")

# Function to plot violin plots
#rm(plots)
rm(violin)
violin <- function(GENE){
		 dat <- plotting_df %>% filter(hugo_ID==GENE)
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
			  annotate(geom = "text", x = 1.5, y = max(plotting_df$value)+1, label=paste0("adj. p value: ", corrected_pval)) +
	   		  ggtitle(paste0("Violin plot of batch1 simulated ", GENE, "\n expression after adjustment via voom()")) +
			  theme(plot.title=element_text(hjust=0.5))
	   		  ggsave(paste0("/scratch/mjpete11/linear_models/linear/simulated_voom_batch/simulated_voom_blocking_with_batch_plots/", GENE, ".png"), device="png")
}
plots <- Map(violin, GENE=SLC)
plots
