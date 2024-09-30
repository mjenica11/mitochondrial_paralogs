# PCA plots of unadjusted, quantile normalize, and combat normalized sweet count data

library(tidyr)
library(data.table)
library(dplyr)
library(stats)
library(ggplot2)
library(ggfortify)
library(tidyverse)
library(stringr)
library(zoo)
library(readxl)
library(janitor)
library(viridis)

# Path
path <- getwd()
path

## Read in filtered not not normalized data 
# Read in clinical data 
nf_counts <- read.table("/scratch/mjpete11/mitochondrial_paralogs/sweet_data/count_matrices/filtered_non_failing_controls.csv", sep=",")
icm_counts <- read.table("/scratch/mjpete11/mitochondrial_paralogs/sweet_data/count_matrices/filtered_ischemic_cardiomyopathy.csv", sep=",")
dcm_counts <- read.table("/scratch/mjpete11/mitochondrial_paralogs/sweet_data/count_matrices/filtered_dilated_cardiomyopathy.csv", sep=",")

# Combine into one dataframe
nf_counts[1:5,1:5]
icm_counts[1:5,1:5]
dcm_counts[1:5,1:5]

# Dimensions
dim(nf_counts) # 24210    16 
dim(icm_counts) # 24688    15  
dim(dcm_counts) # 25571    38 --> 63 samples total

# Make the first row into the colnames
nf_counts <- nf_counts %>% row_to_names(row_number=1) 
nf_counts[1:5,1:5]
colnames(nf_counts)
colnames(nf_counts)[1] <- "index"

icm_counts <- icm_counts %>% row_to_names(row_number=1) 
icm_counts[1:5,1:5]
colnames(icm_counts)
colnames(icm_counts)[1] <- "index"

dcm_counts <- dcm_counts %>% row_to_names(row_number=1) 
dcm_counts[1:5,1:5]
colnames(dcm_counts)
colnames(dcm_counts)[1] <- "index"

# Combine into one dataframe
counts <- list(nf_counts, icm_counts, dcm_counts) %>% reduce(inner_join, by=c("Hugo_ID"))
head(counts)
dim(counts) #  23346 65  --> 63 samples
colnames(counts)
head(counts$Hugo_ID) 

# Drop the index cols
counts$index.x <- NULL
counts$index.y <- NULL
counts$index <- NULL

# Write to file
write.csv(counts, paste0(path, "/combined_filtered_matrix.csv"))

