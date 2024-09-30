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

# path
path <- getwd()
path

#_______________________________________________________________________________ 
# Read in filtered, combined matrix 
#_______________________________________________________________________________ 
# Read in clinical data 
counts <- read.csv("/scratch/mjpete11/mitochondrial_paralogs/sweet_data/count_matrices/combined_filtered_matrix.csv")
head(counts)
dim(counts) # 23346 67 

# Transpose matrix so I can append the covariates to control for
counts[1:5,1:5]
tcounts <- as.data.frame(t(counts))
tcounts[1:5,1:5]
#class(tcounts) # matrix
lapply(counts, class)

genes <- counts$Hugo_ID

# Drop unneccessary rows
tcounts <- counts[-c(1:2),]
tcounts[1:5,1:5]

# Read in sample attributes files
metadata <- read.csv("/scratch/mjpete11/mitochondrial_paralogs/sweet_data/count_matrices/metadata.csv", sep=",")

colnames(metadata)
class(metadata)
dim(metadata) # 64 30

tail(colnames(counts))
head(colnames(counts))
colnames(counts)

# Add a column with the run IDs
meta <- data.frame(run_ID=metadata$Run)
head(meta)
class(meta)
meta <- data.frame(run_ID=meta[!duplicated(meta),])
class(meta)
dim(meta) # 64 1

# Add a columns with the sex in each df
meta$sex <- metadata$"sex" 

# Add a column with the sample ages
meta$age <- metadata$"AGE"

# Add a column with the disease status
meta$disease <- metadata$disease

# Add the sample names
meta$sample_ID <- metadata$BioSample

# Any duplicated rows?
any(duplicated(meta)==TRUE) # FALSE
head(meta)
tail(meta)

# Drop the DCM sample that Salmon could not process
meta <- meta[!meta$sample_ID=="SAMN09484097",]
dim(meta) # 63 5

# Reorder the columns (samples) in the counts df to be in the
# same order as the samples in the organs metadata df (rows)
# Necessary because voom() assumes the rows of the design matrix are the samples
# and they are in the same order as the column of the counts matrix
idx <- match(colnames(counts[,3:ncol(counts)]), meta$sample_ID)
idx
length(idx) # 63
counts <- data.table(counts)
class(counts)
dim(counts) # 23346 65
# Drop index and gene name column
counts$X <- NULL
counts$Hugo_ID <- NULL
counts[1:5,1:5]
ordered_counts <- counts[,..idx]
ordered_counts[1:5,1:5]
class(ordered_counts) # "data.table" "data.frame"
apply(ordered_counts, 2, class) # character
#ordered_counts <- apply(ordered_counts[,3:ncol(counts)], 2, as.numeric)
#apply(ordered_counts, 2, class) # numeric 
dim(ordered_counts) # 23346    61
dim(meta) # 64 5
all(meta$sample_ID==colnames(ordered_counts)) # TRUE

# What is the range
range(ordered_counts) # 0 923219
apply(ordered_counts, 2, range)
head(apply(ordered_counts, 2, range))
apply(ordered_counts, 2, class)

# weird samples: "SAMN09484136","SAMN09484137"

# Make a histogram of the expression distribution after adjusting
png(paste0(path,"/histogram_filtered_counts.png"))
hist(ordered_counts, col="red")
dev.off()

# Check spelling of metadata variables
colnames(meta)
class(meta)

# Check the class of the meta variables
apply(meta, 2, class) # ... everything is a character vector :(

# Change the age class so it is numeric
meta$age <- as.numeric(as.character(meta$age))
class(meta$age) # numeric

# Make design matrix
mat <- model.matrix(~ 0 + meta$sex + meta$age)
nrow(mat)==ncol(ordered_counts) # TRUE
head(mat)
tail(mat)
class(mat)

#_______________________________________________________________________________ 
# Perform quantile normalization 
#_______________________________________________________________________________ 
# Does the number of columns (samples) in the ordered counts matrix
# match the number of rows (samples) in the organs df?
counts_mat <- as.matrix(ordered_counts)
counts_mat <- apply(counts_mat, 2, as.numeric)
apply(ordered_counts, 2, class) # numeric 
head(counts_mat)
ncol(counts_mat)==nrow(meta) # TRUE
nrow(ordered_counts) # 24209 

# Apply voom normalization and quantile normalization
#voom_obj <- voom(counts=counts_mat[1:20,1:20], design=mat[1:20,], normalize.method="quantile", save.plot=TRUE) # TEST 
voom_obj <- voom(counts=counts_mat, design=mat, normalize.method="quantile", save.plot=TRUE) # TEST 

# Adjust for sex and age, but skip quantile normalization
#voom_obj <- voom(counts=counts_mat, design=mat, normalize.method="none", save.plot=TRUE) # TEST 

voom_obj$weights[1:5,1:5] 
ordered_counts[1:5,1:6]

# Make a histogram of the expression distribution after adjusting
png(paste0(path,"/histogram_post_voom_.png"))
#png(paste0(path,"/histogram_post_voom_no_quantile.png"))
hist(voom_obj$weights, col="red")
dev.off()

#dim(voom_obj[["weights"]])==dim(ordered_counts[1:20,1:20]) # TRUE TRUE
dim(voom_obj[["weights"]])==dim(ordered_counts) # TRUE TRUE

# Save the logCPM and voom normalized counts into a separate df
dat <- as.data.frame(voom_obj[["weights"]])
class(dat)

# Add column names back
colnames(dat) <- colnames(voom_obj$E) 
dat[1:5,1:5]

# Add gene name column 
dat$Hugo_ID <- genes
class(dat)

# Move the gene name column to the front
dat <- dat %>% select("Hugo_ID", everything())
dat[1:5,1:5]

# Double check that the column names of the original counts and voom
# transformed counts are in the same order
identical(colnames(dat)[2:ncol(dat)], colnames(counts)) # TRUE

# Write to file
write.csv(dat, "/scratch/mjpete11/mitochondrial_paralogs/sweet_data/qnorm_voom_normalize/voom_qnorm_counts.csv")
#write.csv(dat, "/scratch/mjpete11/mitochondrial_paralogs/sweet_data/qnorm_voom_normalize/voom_only_counts.csv")

# Read to file
voom_df <- read.csv("/scratch/mjpete11/mitochondrial_paralogs/sweet_data/qnorm_voom_normalize/voom_qnorm_counts.csv")
#voom_df <- read.csv("/scratch/mjpete11/mitochondrial_paralogs/sweet_data/qnorm_voom_normalize/voom_only_counts.csv")
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
colnames(meta2)

# Write to file
write.csv(meta2, "/scratch/mjpete11/mitochondrial_paralogs/sweet_data/qnorm_voom_normalize/meta_voom_qnorm_counts.csv")
#write.csv(meta2, "/scratch/mjpete11/mitochondrial_paralogs/sweet_data/qnorm_voom_normalize/meta_voom_only_counts.csv")

