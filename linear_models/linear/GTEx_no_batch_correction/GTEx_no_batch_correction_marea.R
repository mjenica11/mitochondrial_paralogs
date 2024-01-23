# Create dataframes of the GTEx heart and liver data without batch correction

# Libraries
library(plyr)
library(limma)
library(data.table)
library(tidyverse)
library(reshape)
library(stringr)
library(ggplot2)
library(ggpubr)
library(edgeR)

# Read in quantile normalized counts
# Start with 9 samples since I keep getting out of memory errors
counts <- fread("/scratch/mjpete11/mitochondrial_paralogs/linear_models/data/data/filtered_counts.csv", sep=",") # float
counts[1:5,1:5]

# Drop the index column
counts$V1 <- NULL

# Read in sample attributes files
file2 <- read.csv("/scratch/mjpete11/mitochondrial_paralogs/linear_models/data/data/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt", sep="\t")

# Subset heart left ventricle and liver samples into two separate dfs
heart <- file2[file2$SMTSD %in% "Heart - Left Ventricle", ]
liver <- file2[file2$SMTSD %in% "Liver", ]

# Subset the counts by the sample IDs to get dataframes of just the heart and liver counts 
heart2 <- counts %>% select(contains(heart$SAMPID))
liver2 <- counts %>% select(contains(liver$SAMPID))

# Append the hugo IDs back on to the dataframes
heart2$Hugo_ID <- counts$Description
heart2 <- heart2 %>% select(Hugo_ID, everything())
liver2$Hugo_ID <- counts$Description
liver2 <- liver2 %>% select(Hugo_ID, everything())
heart2[1:5,1:5]
liver2[1:5,1:5]

# Write dataframes to file
write.table(heart2, "/scratch/mjpete11/mitochondrial_paralogs/linear_models/data/data/GTEx_heart_filtered_counts.csv",sep=",",row.names=FALSE)
write.table(liver2, "/scratch/mjpete11/mitochondrial_paralogs/linear_models/data/data/GTEx_liver_filtered_counts.csv",sep=",",row.names=FALSE)

