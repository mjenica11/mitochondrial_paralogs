# Perform combat normalization on count data

# Libraries
#install.packages("devtools")
#devtools::install_github("zhangyuqing/sva-devel")
library(data.table)
library(tidyverse)
library(reshape)
library(sva)
library(dplyr)
library(data.table)

# Read in counts
#counts <- fread("/scratch/mjpete11/linear_models/data/combined_simulated_batch1_batch2.csv") # unfiltered
counts_mat <- fread("/scratch/mjpete11/linear_models/data/combined_simulated_001_batch1_batch2.csv") # unfiltered
counts_mat[1:5,1:5]
class(counts_mat)

# Convert to data.frame only and move the ensembl_ID column to the front
counts_mat$V1 <- NULL
#counts_mat$ensembl_ID <- NULL
#counts_mat$hugo_ID <- NULL
counts_mat[1:5,1:5]
class(counts_mat)
counts1 <- as.data.frame(counts_mat)
class(counts1)
dim(counts1) # 61386 2000
counts1[1:5,1:5]
counts1 <- counts1 %>% select(ensembl_ID, everything())

# Moved the filtering step here since I tested different error rate parameters first
# Batch #1 = 0.005
# Batch #2 = 0.001, 0.05, 0.5
# Expression filter: Keep a gene if it has >x counts1 in at least 1 sample 
# Skip the first column because that is just the gene name column
expression_filter <- function(DF, thresh){
						#DAT <- DF[rowSums(DF[,2:ncol(DF)]>=thresh, na.rm=TRUE)>0, ]
						DAT <- DF[rowSums(DF[,2:ncol(DF)]>=thresh, na.rm=TRUE)>0, ]
						return(DAT)
}

# Apply function
counts2 <- expression_filter(DF=counts1, thresh=5)
dim(counts2) # batch #2 error rate = 0.05 --> 60 2001, 0.5 --> 60 2001
counts2[1:5,1:5]

# Round up to integer values
counts2 <- counts2 %>% mutate_if(is.numeric, round)
counts2

# Make design vector for the filtered samples dataframe
batch1_vec <- rep_len('batch1', len=ncol(counts2)/2)
batch2_vec <- rep_len('batch2', len=ncol(counts2)/2)
length(batch1_vec) # 1000
length(batch2_vec) # 1000
batch_vect <- c(batch1_vec, batch2_vec)
length(batch_vect)==ncol(counts2)-1 # TRUE
sapply(batch_vect, class)
batch_vect <- as.factor(batch_vect)
sapply(batch_vect, class)

# Make design matrix with all samples
# since you need a numeric matrix only
#dim(counts_mat) # 53 2000
#batch1_vec <- rep_len('batch1', len=ncol(counts_mat)/2)
#batch2_vec <- rep_len('batch2', len=ncol(counts_mat)/2)
#length(batch1_vec) # 1000
#length(batch2_vec) # 1000
#batch_vect <- c(batch1_vec, batch2_vec)
#mat <- fastDummies::dummy_cols(batch_vect)
#head(mat);tail(mat)
#class(mat) # data.frame
#mat <- data.matrix(mat)
#class(mat) # matrix
#dim(mat) # 2000 3
#nrow(mat);ncol(counts_mat) # 2000 2000 
#nrow(mat)==ncol(counts_mat) # TRUE

# Check if there are any genes left with 0 expression
which(!rowSums(counts2==0))
class(counts2) # data.frame
counts2[1:5,1:5]
counts2[1:5,2:5]
tail(colnames(counts2))

# Combat complaining about offsets with non-finite values
is_finite_data_frame <- function(obj){
		    sapply(obj,FUN = function(x) all(is.finite(x)))
}

all(is_finite_data_frame(counts2[,2:ncol(counts2)]))==TRUE # TRUE
counts2

# Apply ComBat_seq to the data, using parametric empirical Bayesian adjustment
combat_edata <- ComBat_seq(counts=as.matrix(counts2[,2:ncol(counts2)]), batch=batch_vect)
head(combat_edata)
str(combat_edata)

# Append gene names back 
combat <- as.data.frame(combat_edata)
combat[1:5,1:5]
dim(combat)
combat$"ensembl_ID" <- counts2$"ensembl_ID"
head(colnames(combat))
tail(colnames(combat))

# Move the gene names column to the front
combat <- combat %>% select("ensembl_ID", everything())
combat[1:5,1:5]

# Write to file
write.csv(combat, "/scratch/mjpete11/linear_models/data/simulated_error_001_combat_seq.csv")
