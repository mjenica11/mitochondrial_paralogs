# Perform combat normalization on count data

# Libraries
library(data.table)
library(tidyverse)
library(reshape)
library(sva)
library(dplyr)

# Set the path
path <- paste0(getwd(),"/")
path

# Set the count file
file1 <- "non_failing_controls_matrix.csv" 
file2 <- "dilated_cardiomyopathy_matrix.csv"
file3 <- "ischemic_cardiomyopathy_matrix.csv"

# Set the metadata file
metadata <- "metadata.csv"

# Outfile
outfile1 <- "filtered_non_failing_controls.csv"
outfile2 <- "filtered_ischemic_cardiomyopathy.csv"
outfile3 <- "filtered_dilated_cardiomyopathy.csv"

# Read in clinical data 
counts1 <- fread(paste0(path, file1), sep="\t")
counts2 <- fread(paste0(path, file2), sep="\t")
counts3 <- fread(paste0(path, file3), sep="\t")

# Dimensions of the non-failing controls before filtering
# Sample #37 didn't process (the 23rd DCM sample)
# Skip the 37th row (i.e. 23rd DCM sample) when labeling
dim(counts1) # 41219 15 --> 14 non failing  samples
dim(counts2) # 41219 14 --> 13 ischemic samples
dim(counts3) # 41219 37 --> 36 dcm samples
counts1[1:5,1]
counts2[1:5,1:3]
counts3[1:5,1]

# Read in metadata with batch variables
manifest <- read.csv(paste0(path, metadata), header=TRUE, sep = ",")
dim(manifest) # 64 30
head(manifest)

# Samples to drop from metadata bc Salmon couldn't process it: 
# RUN: SRR7426822
# BioSample: SAMN09484097 
manifest <- manifest[!manifest$BioSample=="SAMN09484097",]
dim(manifest) # 63 30

# Change the column names from the run ID to the sample ID
# This is fine because they are uninque IDs
colnames(counts1)[2:ncol(counts1)] <- manifest$BioSample[1:14]
colnames(counts2)[2:ncol(counts2)] <- manifest$BioSample[15:27]
colnames(counts3)[2:ncol(counts3)] <- manifest$BioSample[28:63]

# Change the first column name
colnames(counts1)[1] <- c("Hugo_ID")
colnames(counts2)[1] <- c("Hugo_ID")
colnames(counts3)[1] <- c("Hugo_ID")

dim(counts1)
dim(counts2)
dim(counts3)

# Move the gene names column to the front
counts1 <- counts1 %>% select("Hugo_ID", everything())
counts2 <- counts2 %>% select("Hugo_ID", everything())
counts3 <- counts3 %>% select("Hugo_ID", everything())
counts1[1:5]
counts2[1:5,1]
counts3[1:5,1]

# Expression filter: Keep a gene if it has >x counts in at least 1 sample 
# Skip the first column because that is just the gene name column
expression_filter <- function(DF, thresh){
				DAT <- DF[rowSums(DF[,2:ncol(DF)]>=thresh, na.rm=TRUE)>0, ]
				return(DAT)
}

# Apply function
counts1_2 <- expression_filter(DF=counts1, thresh=5)
dim(counts1_2) # 24209    15 

counts2_2 <- expression_filter(DF=counts2, thresh=5)
dim(counts2_2) #  24687    14

counts3_2 <- expression_filter(DF=counts3, thresh=5)
dim(counts3_2) # 25570    37 

# Use another function to test if it is working...
tmp <- counts2[rowSums(counts2[,2:ncol(counts2)]>=5, na.rm=TRUE)>0,] 
dim(tmp) # 24147 14 

# Write filtered counts to file
write.csv(counts1_2, paste0(path,outfile1))

write.csv(counts2_2, paste0(path,outfile2))

write.csv(counts3_2, paste0(path,outfile3))

# Did any SLCs drop out?
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

# Are any genes missing?
setdiff(SLC, counts1_2$'Hugo_ID')  # non-failing controls: SLC25A2, SLC25A31, and SLC25A47 are missing 
setdiff(SLC, counts2_2$'Hugo_ID')  # ischemic controls: SLC25A2, SLC25A31, and SLC25A47 are missing 
setdiff(SLC, counts3_2$'Hugo_ID')  # dilated controls: SLC25A2, SLC25A31, and SLC25A47 are missing 
