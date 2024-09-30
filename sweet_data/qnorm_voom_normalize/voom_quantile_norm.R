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
library(data.table)
library(janitor)

#_______________________________________________________________________________ 
# Read in combat_seq batch adjusted matrix 
#_______________________________________________________________________________ 
# Read in clinical data 
nf_counts <- read.table("/scratch/mjpete11/mitochondrial_paralogs/sweet_data/count_matrices/filtered_non_failing_controls.csv", sep=",")
icm_counts <- read.table("/scratch/mjpete11/mitochondrial_paralogs/sweet_data/count_matrices/filtered_ischemic_cardiomyopathy.csv", sep=",")
dcm_counts <- read.table("/scratch/mjpete11/mitochondrial_paralogs/sweet_data/count_matrices/filtered_dilated_cardiomyopathy.csv", sep=",")

# Combine into one dataframe
head(nf_counts) 
head(icm_counts) 
head(dcm_counts) 

# Dimensions
dim(nf_counts) # 24148 15
dim(icm_counts) # 25100 24 
dim(dcm_counts) # 25100 31 

# Combine into one dataframe
counts <- list(nf_counts, icm_counts, dcm_counts) %>% reduce(inner_join, by="V1")
head(counts)
dim(counts) # 24148 68 

# Make the first row into the colnames
counts <- counts %>% row_to_names(row_number=1) 
dim(counts) # 24147 68 
head(counts)

# Transpose matrix so I can append the covariates to control for
counts[1:5,1:5]
#tcounts <- t(counts)
#tcounts[1:5,1:5]
#class(tcounts) # matrix
lapply(counts, class)

# Read in sample attributes files
file <- read.csv("/scratch/mjpete11/mitochondrial_paralogs/sweet_data/count_matrices/metadata.csv", sep=",")

colnames(file)
class(file)
dim(file) # 64 30

tail(colnames(counts))
head(colnames(counts))
colnames(counts)

# drop the irrelevant column names
colnames <- colnames(counts)[-c(1,2,16,39)]
colnames

# Add a column with the run IDs
meta <- data.frame(run_ID=colnames)
head(meta)
class(meta)
meta <- data.frame(run_ID=meta[!duplicated(meta),])
class(meta)

# Add a columns with the sex in each df
meta$sex <- file$"sex" 

# Add a column with the biosample IDs 
meta$run_ID <- file$"Run"
class(meta$run_ID)

# Add a column with the sample ages
meta$age <- file$"AGE"

# Add a column with the disease status
meta$disease <- file$disease

# Add the sample names
meta$sample_ID <- file$BioSample

# Any duplicated rows?
any(duplicated(meta)==TRUE) # FALSE
head(meta)
tail(meta)

# Drop the irrelevant columns
counts <- counts[,-c(1,2,16,39)]
colnames(counts)
class(counts)

# Reorder the columns (samples) in the counts df to be in the
# same order as the samples in the organs metadata df (rows)
# Necessary because voom() assumes the rows of the design matrix are the samples
# and they are in the same order as the column of the counts matrix
idx <- match(meta$sample_ID, colnames(counts))
counts <- as.data.table(counts)
ordered_counts <- counts[,..idx]
all(meta$sample_ID==colnames(ordered_counts)) # TRUE

# Make design matrix
mat <- model.matrix(~ 0 + meta$disease)
nrow(mat)==ncol(ordered_counts) # TRUE

#_______________________________________________________________________________ 
# Perform quantile normalization 
#_______________________________________________________________________________ 
# Does the number of columns (samples) in the ordered counts matrix
# match the number of rows (samples) in the organs df?
counts_mat <- as.matrix(ordered_counts)
counts_mat <- apply(counts_mat, 2, as.numeric)
head(counts_mat)
ncol(counts_mat)==nrow(meta) # TRUE
nrow(ordered_counts) # 24147 

print("line 170")
# Apply voom normalization and quantile normalization
#voom_obj <- voom(counts=counts_mat[1:20,1:20], design=mat[1:20,], normalize.method="quantile", save.plot=TRUE) # TEST 
voom_obj <- voom(counts=counts_mat, design=mat, normalize.method="quantile", save.plot=TRUE) # TEST 

voom_obj$weights[1:5,1:5] 
ordered_counts[1:5,1:6]

#dim(voom_obj[["weights"]])==dim(ordered_counts[1:20,1:20]) # TRUE TRUE
dim(voom_obj[["weights"]])==dim(ordered_counts) # TRUE TRUE

print("line 181")
# Save the logCPM and voom normalized counts into a separate df
dat <- as.data.frame(voom_obj[["weights"]])

print("line 185")
# Add column names back
colnames(dat) <- colnames(voom_obj$E) 

print("line 189")
# Add gene name column
#gene_names <- organs$"gene"[1:20]
gene_names <- nf_counts$V2 
gene_names

# Add gene name column 
ordered_counts$Hugo_ID <- gene_names[2:length(gene_names)]

# Move the gene name column to the front
ordered_counts <- ordered_counts %>% select("Hugo_ID", everything())

ordered_counts[1:5,1:5]

# Double check that the column names of the original counts and voom
# transformed counts are in the same order
identical(colnames(ordered_counts)[2:ncol(ordered_counts)], colnames(counts)) # TRUE

# Write to file
write.csv(ordered_counts, "/scratch/mjpete11/mitochondrial_paralogs/sweet_data/qnorm_voom_normalize/voom_qnorm_counts.csv")

# Read to file
voom_df <- read.csv("/scratch/mjpete11/mitochondrial_paralogs/sweet_data/qnorm_voom_normalize/voom_qnorm_counts.csv")
names(voom_df) <- gsub(x=names(voom_df), pattern="\\.", replacement="-")
head(voom_df)
voom_df[1:3,1:3]

# Convert column names (sample ID) into column 
voom_df$X <- NULL
voom_df2 <- melt(voom_df)
head(voom_df2)
dim(voom_df2)
voom_df2[1:3,1:3]
colnames(voom_df2)[2] <- "sample_ID"
colnames(voom_df2)[3] <- "log2_cpm"

# Add voom quantile normalized counts to organs df for easy plotting
meta2 <- left_join(meta, voom_df2, by=c("sample_ID"))
meta2[1:5,1:7]

# Write to file
write.csv(meta2, "/scratch/mjpete11/mitochondrial_paralogs/sweet_data/qnorm_voom_normalize/voom_qnorm_counts.csv")
