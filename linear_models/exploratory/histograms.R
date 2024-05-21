# Purpose: Linear models to asses if SLC paralogs are dependent on energetic state?

# Libraries
library(pscl)
library(MASS)
library(plyr)
library(dplyr)
library(data.table)
library(tidyverse)
library(reshape)
library(stringr)
library(ggplot2)
library(ggpubr)
library(ggtext)
library(rstatix)
library(stats)
library(scales)
library(edgeR)

# Read in combat_seq normalized counts 
#counts <- fread("/scratch/mjpete11/mitochondrial_paralogs/linear_models/data/data/combat_seq_filtered.csv", sep=",")
counts <- fread("/scratch/mjpete11/mitochondrial_paralogs/linear_models/data/data/voom_quantile_normalized_counts.csv", sep=",")
counts[1:5,1:5]

# Drop the index column
counts$V1 <- NULL

# Generate the design matrix
# Samples are rows and columns are covariates of interest (heart and liver)

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

# List of metabolic biomarkers 
biomarkers <- c("SLC16A1", "SLC16A2", "HK1", "HK2", "HK3", "GPI", "PFKM","PFKX", "PFKP", "ALDOA", "ALDOB", "TPI1", "GAPDH", "PGK1","PGK2", "PGAM1", "PGAM5", "PGAM2", "PGAM4", "ENO1","PKLR", "LDHA", "LDHB", "LDHC", "LDHD", "CACT", "CPT1A","CPT1B", "CPT2", "ACADVL", "ACADL", "ACADM", "ACADS","ACADSB", "ACAD11", "ACAD8", "ACAD10", "ACAD9", "ECHS1","ECH1", "ECI1", "ECI2", "ECHS1", "HADHA", "PDHA1", "PDHA2","PDHB", "PDP1", "PDHX", "CS", "ACO1", "ACO2", "IDH3G","IDH1", "IDH2", "IDH3A", "IDH3B", "OGDH", "SUCLG2","SUCLG1", "SUCLA2", "SDHB", "SDHD", "SDHA", "SDHC", "SDHAF2","FH", "MDH1", "MDH2", "MT-ND1")

genes <- c(SLC, biomarkers)
genes

# Subset the SLC25 genes
counts <- as.data.frame(counts)
sub_df <- counts[counts$"Description" %in% genes, ]

# Write to file
#write.table(sub_df, "/scratch/mjpete11/linear_models/data/SLC_df_voom_combat_seq.csv", sep=",")

# Read in count df of just SLC genes
#sub_df <- read.csv("/scratch/mjpete11/mitochondrial_paralogs/linear_models/data/data/SLC_df_voom_combat_seq.csv", sep=",")
#names(sub_df) <- gsub(x=names(sub_df), pattern="\\.", replacement="-")

# All present!
#setdiff(SLC, sub_df$'Description') 
setdiff(genes, sub_df$'Description')  # PFKX PGAM2 CACT

# Read in sample attributes files
file <- read.csv("/scratch/mjpete11/mitochondrial_paralogs/linear_models/data/data/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt", sep="\t")

# Subset heart left ventricle and liver samples into two separate dfs
heart <- file[file$SMTSD %in% "Heart - Left Ventricle", ]
liver <- file[file$SMTSD %in% "Liver", ]

# Subset the SLC25 gene count df by the sample IDs that match the IDs in the 
# GTEx sample annotation df
heart2 <- sub_df %>% select(contains(heart$SAMPID))
liver2 <- sub_df %>% select(contains(liver$SAMPID))

# Append the gene name
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

# Any duplicated rows?
any(duplicated(heart4)==TRUE) # FALSE
any(duplicated(liver4)==TRUE) # FALSE

# Change the name of the 'variable' column to 'SAMPID' to match columns
colnames(heart4)[3] <- "SAMPID" 
colnames(liver4)[3] <- "SAMPID" 

# Add a subject ID column
heart4$SUBJID <- str_sub(heart4$SAMPID, start=1L, end=10L)
liver4$SUBJID <- str_sub(liver4$SAMPID, start=1L, end=10L)

# Add column with the expression batch ID (SMGEBTCH) and the type of genotype
# or expression batch (SMGEBTCHT), the RIN number (SMRIN) and tissue type (SMTSISCH)
heart4$SMGEBTCH <- file$SMGEBTCH[match(heart4$SAMPID, file$SAMPID)]
liver4$SMGEBTCH <- file$SMGEBTCH[match(liver4$SAMPID, file$SAMPID)]

heart4$SMGEBTCHT <- file$SMGEBTCHT[match(heart4$SAMPID, file$SAMPID)]
liver4$SMGEBTCHT <- file$SMGEBTCHT[match(liver4$SAMPID, file$SAMPID)]

heart4$SMRIN <- file$SMRIN[match(heart4$SAMPID, file$SAMPID)]
liver4$SMRIN <- file$SMRIN[match(liver4$SAMPID, file$SAMPID)]

heart4$SMTSISCH <- file$SMTSISCH[match(heart4$SAMPID, file$SAMPID)]
liver4$SMTSISCH <- file$SMTSISCH[match(liver4$SAMPID, file$SAMPID)]

# Keep heart and liver samples from all individuals
organs <- rbind(heart4, liver4)
length(unique(organs$SUBJID)) # 508 individuals

# Are any rows duplicated?
any(duplicated(organs)) # FALSE

# Change the name of 'value' column to 'combat_seq_counts'
#colnames(organs)[4] <- "combat_seq_counts"
colnames(organs)[4] <- "voom_qnorm_counts"

# Convert combat_seq adjusted counts to log2 and CPM
# Values of zero will be set to 2
#organs$log2_cpm <- cpm(organs$combat_seq_counts, log=TRUE, prior.count=2)
organs$log2_cpm <- cpm(organs$voom_qnorm_counts, log=TRUE, prior.count=2)

# Write organs df to file
#write.table(organs, "/scratch/mjpete11/linear_models/data/organs_combat_seq_all_heart_liver_samples.csv", sep=",")

# Set y lim
#range(organs$combat_seq_counts) # c(0 265290) 
range(organs$log2_cpm) # c(-5.96, 11.06) 

# Why does SLC25A47 look like it only has zero values? 
# Make layered histograms clearer so I can see the other layer better
A47 <- organs[organs$gene=="SLC25A47",]
range(A47$combat_seq_counts) # 0 67810

# Histogram plot function
rm(plots)
rm(histogram_plots)
# Negative binomial 2SRI violin and jitter plot function
histogram_plots <- function(GENE){
	dat <- organs %>% filter(gene == GENE)
    # Violin plot
	#p <- ggplot(dat, aes(x = combat_seq_counts, fill=organ)) +
	p <- ggplot(dat, aes(x = voom_qnorm_counts, fill=organ)) +
	geom_histogram(bins=100) +
	scale_color_manual(values = c("green", "purple")) +
#	geom_bin2d(bins = 100) +
	#ylim(c(0,30)) +
	ylab("frequency") +
#	xlab("filtered and combat_seq adjusted counts") +
	xlab("filtered, quantile adjusted, and voom adjusted counts") +
	theme(plot.title = element_textbox_simple(size=10)) + # automatically wrap long titles
	ggtitle(paste0("Histogram plot of ", GENE, " expression in heart and liver")) 
	ggsave(paste0("/scratch/mjpete11/mitochondrial_paralogs/linear_models/data/data/histogram_plots/", GENE, ".png"), p, device="png")
}
plot.new()
plots <- Map(histogram_plots, GENE=genes)
dev.off()
