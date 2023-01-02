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
# Read in combat_seq batch adjusted matrix 
#_______________________________________________________________________________ 
# Read in GTEx counts
counts <- fread("/scratch/mjpete11/linear_models/data/combat_seq_filtered.csv", sep=",")

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
sub_df <- read.csv("/scratch/mjpete11/linear_models/data/SLC_df_voom_combat_seq.csv", sep=",")
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
colnames(tmp1) <- c("SUBJID", "gene", "organ", "SAMPID", "value", "SMRIN", "SMTSISCH", "SMGEBTCH", "SMGEBTCHT")
colnames(tmp2) <- c("SUBJID", "gene", "organ", "SAMPID", "value", "SMRIN", "SMTSISCH", "SMGEBTCH", "SMGEBTCHT")
organs0 <- rbind(tmp1, tmp2)

# Are any rows duplicated?
duplicated(organs0) # Yes

# Drop duplicated rows
organs <- organs0[!duplicated(organs0),]

# Make each sample ID unique
# There will be 53 gene expression measurements from each sample ID, so
# add increasing intergers to make a unique gene x sample ID
organs$SAMPID <- make.unique(organs$SAMPID, sep="_")

# Write organs df to file
write.table(organs, "/scratch/mjpete11/linear_models/data/organs_combat_seq_voom.csv", sep=",")

# Use this file to create the design matrix
organs <- read.csv("/scratch/mjpete11/linear_models/data/organs_combat_seq_voom.csv", sep=",")

# Drop samples from the SLC counts df that are not from paired heart or liver
keeps <- as.vector(organs$SAMPID)
SLC_heart_liver_counts <- subset(sub_df, select=keeps) 

# Are there the same number of samples in the counts df and the organs metadata df? 
ncol(SLC_heart_liver_counts)==nrow(organs) # TRUE

# Drop the decimals introduced into the column names
#colnames(SLC_heart_liver_counts) <- keeps

# Reorder the columns (samples) in the SLC_heart_liver_counts df to be in the
# same order as the samples in the organs metadata df (rows)
# Necessary because voom() assumes the rows of the design matrix are the samples
# and they are in the same order as the column of the counts matrix
idx <- match(organs$SAMPID, colnames(SLC_heart_liver_counts))
SLC_heart_liver_counts <- as.data.table(SLC_heart_liver_counts)
ordered_counts <- SLC_heart_liver_counts[,..idx]
all(organs$SAMPID==colnames(ordered_counts)) # TRUE

# Make design matrix
mat <- model.matrix(~ 0 +organs$organ)
nrow(mat)==ncol(ordered_counts) # TRUE

#_______________________________________________________________________________ 
# Perform quantile normalization 
#_______________________________________________________________________________ 
# Does the number of columns (samples) in the ordered counts matrix
# match the number of rows (samples) in the organs df?
counts_mat <- as.matrix(ordered_counts)
ncol(counts_mat)==nrow(organs) # TRUE
nrow(ordered_counts)

# Apply voom normalization and quantile normalization
#voom_obj <- voom(counts=counts_mat[1:20,1:20], design=mat[1:20,], normalize.method="quantile", save.plot=TRUE) # TEST 
voom_obj <- voom(counts=counts_mat, design=mat, normalize.method="quantile", save.plot=TRUE) # TEST 

voom_obj$weights[1:5,1:5] 
ordered_counts[1:5,1:6]

#dim(voom_obj[["weights"]])==dim(ordered_counts[1:20,1:20]) # TRUE TRUE
dim(voom_obj[["weights"]])==dim(ordered_counts) # TRUE TRUE

# Save the logCPM and voom normalized counts into a separate df
dat <- as.data.frame(voom_obj[["weights"]])

# Add column names back
colnames(dat) <- colnames(voom_obj$E) 

# Add gene name column
#gene_names <- organs$"gene"[1:20]
gene_names <- sub_df$'Description' 

# Make column names unique since there are 53 gene expression measurements
# from the same sample ID
samples <- make.unique(colnames(dat), sep="_")
colnames(dat) <- samples

# Add gene name column 
dat$"Description" <- gene_names

# Move the gene name column to the front
dat2 <- dat %>% select("Description", everything())

dat2[1:5,1:5]
ordered_counts[1:5,1:5]

# Double check that the column names of the original counts and voom
# transformed counts are in the same order
identical(colnames(ordered_counts), colnames(dat2[,2:ncol(dat2)])) # TRUE

# Write to file
write.csv(dat2, "/scratch/mjpete11/linear_models/data/voom_quantile_normalized_counts2.csv")

# Read to file
voom_df <- read.csv("/scratch/mjpete11/linear_models/data/voom_quantile_normalized_counts2.csv")
names(voom_df) <- gsub(x=names(voom_df), pattern="\\.", replacement="-")

# Convert column names (sample ID) into column 
voom_df$X <- NULL
voom_df2 <- melt(voom_df)
colnames(voom_df2)[2] <- "SAMPID"
colnames(voom_df2)[1] <- "gene"

# Add voom quantile normalized counts to organs df for easy plotting
organs2 <- left_join(organs, voom_df2, by=c("gene", "SAMPID"))
colnames(organs2)[10] <- "voom"
colnames(organs2)[5] <- "counts"

# Write to file
write.csv(organs2, "/scratch/mjpete11/linear_models/data/voom_quantile_normalized_counts3.csv")
