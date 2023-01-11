# Try quantile normalization and logCPM transformation with limma voom

# Libraries
library(data.table)
library(tidyverse)
library(reshape)
library(limma)
library(plyr)
library(stringr)
library(ggplot2)
library(ggpubr)

#_______________________________________________________________________________ 
# Read in filtered counts 
#_______________________________________________________________________________ 
counts <- fread("/scratch/mjpete11/linear_models/data/filtered_counts.csv", sep=",")

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

# Subset the SLC25 genes
counts <- as.data.frame(counts)
sub_df <- counts[counts$"Description" %in% SLC, ]

# Write to file
#write.table(sub_df, "/scratch/mjpete11/linear_models/data/SLC_df_voom_combat_seq.csv", sep=",")

# Read in count df of just SLC genes
#sub_df <- read.csv("/scratch/mjpete11/linear_models/data/SLC_df_voom_combat_seq.csv", sep=",")
names(sub_df) <- gsub(x=names(sub_df), pattern="\\.", replacement="-")

# All present!
setdiff(SLC, sub_df$'Description') 

# Read in sample attributes files
file <- read.csv("/scratch/mjpete11/linear_models/data/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt", sep="\t")

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

# Dataframe with only samples from people who donated both a heart and liver: organs
# e.g. paired samples
organs_tmp <- merge(heart4, liver4, by="SUBJID") 
length(unique(organs_tmp$SUBJID)) # 148 individuals
tmp1 <- organs_tmp[,c("SUBJID", "gene.x", "organ.x", "SAMPID.x", "value.x", "SMRIN.x", "SMTSISCH.x", "SMGEBTCH.x", "SMGEBTCHT.x")]
tmp2 <- organs_tmp[,c("SUBJID", "gene.y", "organ.y", "SAMPID.y", "value.y", "SMRIN.y", "SMTSISCH.y", "SMGEBTCH.x", "SMGEBTCHT.x")]
colnames(tmp1) <- c("SUBJID", "gene", "organ", "SAMPID", "filtered_counts", "SMRIN", "SMTSISCH", "SMGEBTCH", "SMGEBTCHT")
colnames(tmp2) <- c("SUBJID", "gene", "organ", "SAMPID", "filtered_counts", "SMRIN", "SMTSISCH", "SMGEBTCH", "SMGEBTCHT")
organs0 <- rbind(tmp1, tmp2)

# Are any rows duplicated?
duplicated(organs0) # Yes

# Drop duplicated rows
organs <- organs0[!duplicated(organs0),]
length(unique(organs_tmp$SUBJID)) # 148 individuals

# Write organs df to file
write.table(organs, "/scratch/mjpete11/linear_models/data/filtered_counts_organ_metadata.csv", sep=",")

