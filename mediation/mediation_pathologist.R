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
tail(act)
class(act) # matrix array
act <- as.data.frame(act)
class(act) # data.frame
lapply(act, class)
dat <- mutate_all(act, function(x) as.numeric(x))
all(lapply(dat, class)=="numeric") # TRUE
head(dat)
tail(dat)
class(dat)
dim(dat) # 296 292
any(as.logical((sapply(dat, is.na))))==TRUE # FALSE 

# Change the rownames to be the first column
dat$SAMPID <- rownames(dat)

# Move the last column to the first
dat <- dat %>% select(SAMPID, everything())
dat[1:5,1:5]

# Drop the rownames column
rownames(dat) <- NULL
dat[1:5,1:5]

# Add the dativity score per sample to the organs df by merging the dataframes
dat1 <- merge(dat, organs, by="SAMPID")
dat1[1:5,1:5]
tail(dat1)
dim(dat1) # 4191656 
dim(dat) #296 292
dim(organs) # 4191656 7

# Drop the duplicate rows
any(duplicated(dat1)) # TRUE
dat2 <- dat1[!duplicated(dat1),]
any(duplicated(dat2)) # FALSE 
dim(dat2) # 35224 298

# Mediation analysis
# Every pathway ~ every gene + rna integrity number + total ischemic time

# Step 1: Estimate the total effect
act1$

# Step 2: Path A (X on M)

# Step 3:Path B (M on Y, controlling for X)

# Step 4: Reversed Path C (Y on X, controlling for M)

# Summary

