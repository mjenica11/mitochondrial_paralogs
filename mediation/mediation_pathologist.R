# Script to do mediation analysis on the activity scores 

# Library
library(data.table)
library(janitor)
library(dplyr)

# Read in the pathways and activity scores at 50% activity
activity <- fread("/scratch/mjpete11/mitochondrial_paralogs/linear_models/data/data/filtered_activity_scores_50percent.csv",sep=",")
activity[1:5,1:5]

# Read in the gene counts with the RIN and activity scores
organs <- read.csv("/scratch/mjpete11/mitochondrial_paralogs/linear_models/data/data/organs_biomarkers_combat_seq3.csv", sep=",")
organs[1:5,1:5]

# Change the .s to  _ in the sample names in the activity dataframe
colnames(activity) <- gsub(x=colnames(activity), pattern="\\.", replacement="-")

# Transpose the activity df 
activity <- t(activity)

# Move the first row to the column name
act <- activity %>% row_to_names(row_number=1)
act[1:5,1:5]
class(act) # matrix array

# Convert the matrix to a dataframe
act <- as.data.frame(act)
class(act) # data.frame 

# Change the rownames to be the first column
act$SAMPID <- rownames(act)

# Move the last column to the first
act <- act %>% select(SAMPID, everything())

# Drop the rownames column
rownames(act) <- NULL

# Add the activity score per sample to the organs df by merging the dataframes
act1 <- merge(act, organs, by="SAMPID")
act1[1:5,1:5]

# Mediation analysis
# Every pathway ~ every gene + rna integrity number + total ischemic time

# Step 1: Estimate the total effect


# Step 2: Path A (X on M)

# Step 3:Path B (M on Y, controlling for X)

# Step 4: Reversed Path C (Y on X, controlling for M)

# Summary

