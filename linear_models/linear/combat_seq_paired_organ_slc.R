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

# Read in quantile normalized counts
# Start with 9 samples since I keep getting out of memory errors
counts <- fread("/scratch/mjpete11/linear_models/data/quantile_normalized_combat_seq.csv", sep=",") # float

# Drop the index column
counts$V1 <- NULL

# List of SLC25 paralogs
SLC <- c("SLC25A1", "SLC25A2", "SLC25A3", "SLC25A4", "SLC25A5", "SLC25A6",
		 "SLC5A7", "UCP2", "UCP3", "SLC25A10", "SLC25A11", "SLC25A12", 
		 "SLC25A13", "SLC25A14", "SLC25A15", "SLC25A16", "SLC25A17", 
		 "SLC25A18", "SLC25A19", "SLC25A20", "SLC25A21", "SLC25A22", 
		 "SLC25A23", "SLC25A24", "SLC25A25", "SLC25A26", "SLC25A27", 
		 "SLC25A28", "SLC25A29", "SLC25A30", "SLC25A3", "SLC25A32", 
		 "SLC25A33", "SLC25A34", "SLC25A35", "SLC25A36", "SLC25A37", 
		 "SLC25A38", "SLC25A39", "SLC25A40", "SLC25A41", "SLC25A42", 
		 "SLC25A43", "SLC25A44", "SLC25A45", "SLC25A46", "SLC25A47",
		 "SLC25A48", "MTCH1", "MTCH2", "SLC25A51", "SLC25A52", "SLC25A53")

# Subset the SLC25 genes
counts <- as.data.frame(counts)
sub_df <- counts[counts$"Description" %in% SLC, ]

# Write to file
write.table(sub_df, "/scratch/mjpete11/linear_models/data/sub_df_combat_seq.csv", sep=",")

# A2, A4, A12, A26, A33, A34, A38, A44, A46, A48 were filtered out
# because they had <5 counts across all samples
setdiff(SLC, sub_df$'Description') 

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
organs_tmp <- merge(heart4, liver4, by="SUBJID") 
length(unique(organs_tmp$SUBJID)) # 148 individuals
tmp1 <- organs_tmp[,c("SUBJID", "gene.x", "organ.x", "SAMPID.x", "value.x", "SMRIN.x", "SMTSISCH.x")]
tmp2 <- organs_tmp[,c("SUBJID", "gene.y", "organ.y", "SAMPID.y", "value.y", "SMRIN.y", "SMTSISCH.y")]
colnames(tmp1) <- c("SUBJID", "gene", "organ", "SUBJID", "value", "SMRIN", "SMTSISCH")
colnames(tmp2) <- c("SUBJID", "gene", "organ", "SUBJID", "value", "SMRIN", "SMTSISCH")
organs <- rbind(tmp1, tmp2)
colnames(organs)[4] <- "SAMPID"

# Write organs df to file
#write.table(organs, "/scratch/mjpete11/linear_models/data/organs_combat_seq.csv", sep=",")

# Read organs df back in
#organs <- read.csv("/scratch/mjpete11/linear_models/data/organs.csv", sep=",")

# Histogram to visualize distribution of values after quantile normalization
###### TEST ######
# Subset one gene from the organ df
A1 <- organs %>% filter(gene == "SLC25A1")
dA1 <- distinct(na.omit(A1))
dA1$value <- round(dA1$value, 2) 
p <- ggplot(dA1, aes(value)) +
     geom_histogram(bins=100) +
     ylab("frequency") +
     xlab("counts") +
     ggtitle(paste0("Histograms of quantile normlized ", "SLC25A1", " expression")) 
ggsave("test.pdf", p, device="pdf")
###### TEST ######

# Historgam function
histogram <- function(GENE){
     dat <- organs %>% filter(gene == GENE)
     dat2 <- distinct(na.omit(dat))
	 dat2$value <- round(dat2$value, 2)
     p <- ggplot(dat2, aes(value)) +
       geom_histogram(bins=100) +
       ylab("frequency") +
       xlab("counts") +
       ggtitle(paste0("Histograms of combat_seq quantile normlized ", GENE, " expression")) 
	 ggsave(paste0("/scratch/mjpete11/linear_models/linear/combat_seq_histograms/", GENE, ".pdf"), p, device="pdf")
#     return(p)
}
plts <- Map(histogram, GENE=SLC)

# Regress out effect of RIN and ischemic time and remove from the model 
# before making the violin plots

###### TEST ######
# Double check that for each value in the subject ID column there is a 
# corresponding heart and liver value in the organ column
any(tapply(dA1$organ, dA1$SUBJID, function(x) length(unique(x))) != 2) #FALSE
# Are the number of rows divided by 2 equal to the number of unique subject IDs?
nrow(dA1)/2==148 # TRUE
# Regress out the effect due to RIN and ischemic time
model1 <- lm(value ~ SMRIN + SMTSISCH, data=dA1)
# Fit model
dA1$resids <- residuals(model1)
# Linear regression on organ
model2 <- lm(value ~ organ + resids, data=dA1)
# Add fitted SLC values to df for plotting
dA1$fitted_values <- model2$fitted.values

# Test violin jitter plot function
p <- ggplot(dA1, aes(x = organ, y = fitted_values)) +
	 geom_violin(aes(colour = organ)) +
	 geom_jitter(aes(colour = organ))  +
	 ylab("fitted value") +
	 xlab("organ") +
	 ggtitle(paste0("Jitter plot of SLC25A1 expression between heart and liver after 2SRI")) +
#   scale_y_continuous(breaks = round(seq(min(dA1$fitted_values), max(dA1$fitted_values), by = 2000),0)) +
	 stat_compare_means(method = "t.test") +
	 stat_summary(fun.data = "mean_cl_boot", geom = "crossbar", colour = "red", width = 0.2)
p
ggsave(paste0("SLC25A1_violin_plot.pdf"))
###### TEST ######

SLCs <- c("SLC25A1", "SLC25A3", "SLC25A5", "SLC25A6",
		 "SLC5A7", "UCP2", "UCP3", "SLC25A10", "SLC25A11",  
		 "SLC25A13", "SLC25A14", "SLC25A15", "SLC25A16", "SLC25A17", 
		 "SLC25A18", "SLC25A19", "SLC25A20", "SLC25A21", "SLC25A22", 
		 "SLC25A23", "SLC25A24", "SLC25A25", "SLC25A27", 
		 "SLC25A28", "SLC25A29", "SLC25A30", "SLC25A3", "SLC25A32", 
		 "SLC25A35", "SLC25A36", "SLC25A37", 
		 "SLC25A39", "SLC25A40", "SLC25A41", "SLC25A42", 
		 "SLC25A43", "SLC25A45", "SLC25A47", "MTCH1", 
		 "MTCH2", "SLC25A51", "SLC25A52", "SLC25A53")

# Violin and jitter plot function
violin_plots <- function(GENE){
			dat <- organs %>% filter(gene == GENE)
			dat2 <- na.omit(dat)
			dat3 <- distinct(na.omit(dat2))
			# Regress out effect of RIN and ischemic time on SLC
			model1 <- lm(value ~ SMRIN + SMTSISCH, data=dat3)
			# Add fitted values to dataframe
			dat3$resids <- residuals(model1)
			# Regress organ on SLC using previous fitted regression values as offset
			model2 <- lm(value ~ organ + resids, data=dat3)
			# Add residuals from second linear model to dataframs
			dat3$fitted_values <- fitted.values(model2)
			# Violin plot
			p <- ggplot(dat3, aes(x = organ, y = value)) +
			geom_violin(aes(colour = organ)) +
			geom_jitter(aes(colour = organ))  +
            # scale_y_continuous(breaks = round(seq(min(dat3$fitted_values), max(dat3$fitted_values), by = 5000),0)) +
			ylab("counts") +
			xlab("organ") +
			ggtitle(paste0("Violin plot of ", GENE, " expression between heart and liver after 2SRI")) +
			stat_compare_means(method="t.test") +
			stat_summary(fun.data = "mean_cl_boot", geom = "crossbar", colour = "red", width = 0.2)
	        ggsave(paste0("/scratch/mjpete11/linear_models/linear/combat_seq_paired_violin_plots/", GENE, ".pdf"), p, device="pdf")
}
plots <- Map(violin_plots, GENE=SLCs)
plots[[1]]
plots

