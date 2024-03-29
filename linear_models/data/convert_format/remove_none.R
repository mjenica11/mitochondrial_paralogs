# Figure out how many metabolic reactions are encoded as "none"
# then drop them from the heart and liver RAS (reaction activity scores)

# Load libraries
library(data.table)
library(dplyr)
library(stringr)

# Read in heart and liver RAS scores generated by Expression2Ras galaxy server
heart_RAS <- fread("/scratch/mjpete11/mitochondrial_paralogs/linear_models/data/data/heart2_RAS_recon22.txt")
class(heart_RAS) # data.table data.frame
heart_RAS <- as.data.frame(heart_RAS)
class(heart_RAS) # data.frame
heart_RAS[1:5,1:5]
# Check how did fread() read in the cells encoded as "none"
heart_RAS["RE1240C",1:5] # <NA>

# How many metabolic reactions are encoded as "none"?
heart_RAS
class(heart_RAS) # data.frame
heart_RAS <- as.data.frame(heart_RAS)
class(heart_RAS) # data.frame
heart_RAS["RE1240C",1:5] # <NA>
any(is.na(heart_RAS))==TRUE # It doesn't pick up on NA values if the rowname is NA
# Drop rows with NA as rhe rowname
heart2 <- heart_RAS[rowSums(heart_RAS == "None")==0, , drop = FALSE]
heart2["RE1240C",1:5] # <NA>
heart2[8:10,1:5] # <NA>
any(is.na(heart2))==TRUE # It doesn't pick up on NA values if the rowname is NA

# Repeat with liver
liver_RAS <- fread("/scratch/mjpete11/mitochondrial_paralogs/linear_models/data/data/liver_RAS_Recon22.txt")
class(liver_RAS) # data.table data.frame
liver_RAS <- as.data.frame(liver_RAS)
class(liver_RAS) # data.frame
liver_RAS[1:5,1:5]
# Check how did fread() read in the cells encoded as "none"
liver_RAS["RE1240C",1:5] # <NA>

# How many metabolic reactions are encoded as "none"?
liver_RAS
class(liver_RAS) # data.frame
liver_RAS <- as.data.frame(liver_RAS)
class(liver_RAS) # data.frame
liver_RAS["RE1240C",1:5] # <NA>
any(is.na(liver_RAS))==TRUE # It doesn't pick up on NA values if the rowname is NA
# Drop rows with NA as rhe rowname
liver2 <- liver_RAS[rowSums(liver_RAS == "None")==0, , drop = FALSE]
liver2["RE1240C",1:5] # <NA>
liver2[8:10,1:5] # <NA>
any(is.na(liver2))==TRUE # It doesn't pick up on NA values if the rowname is NA

# Check that equal numbers of rows were dropped in the heart and liver RAS dataframes
dim(heart2)==dim(liver2) # TRUE TRUE
dim(heart2) # 4741 149
dim(liver2) # 4741 149

# Sort the columns to be in the same order
# The sample names will only be partial matches because the last 5 characters indicate
# organ that was sampled

# Create a vector of strings of the sample IDs from the column names
heart_sample_id <- str_sub(colnames(heart2), start=1L, end=10L)
liver_sample_id <- str_sub(colnames(liver2), start=1L, end=10L)

# Sort the vector of sample IDs to be in the same order
heart_sampleID_sorted <- heart_sample_id[order(heart_sample_id)]
liver_sampleID_sorted <- liver_sample_id[order(heart_sample_id)]

# Check the vectors are identical
all.equal(heart_sampleID_sorted, liver_sampleID_sorted) # TRUE

# Move the "reactions" column to the front
heart_colnames <- c(heart_sampleID_sorted[149], heart_sampleID_sorted)
heart_colnames <- heart_colnames[-length(heart_colnames)]
length(heart_colnames)==ncol(heart2) # TRUE

liver_colnames <- c(liver_sampleID_sorted[149], liver_sampleID_sorted)
liver_colnames <- liver_colnames[-length(liver_colnames)]
length(liver_colnames)==ncol(liver2) # TRUE

# Check the vectors are identical
all.equal(heart_colnames, liver_colnames) # TRUE

# Replace the original column names with the sample IDs only
# I am guessing the samples need to be in the same order for MaREA
class(heart2) # data.frame
setDT(heart2)
class(heart2) # data.table data.frame
idx <- match(heart_colnames, liver_colnames)
setcolorder(heart2, match(seq_along(heart2), idx))
setcolorder(liver2, match(seq_along(liver2), idx))
head(colnames(liver2));head(colnames(heart2)) # Sample IDs are in the right order
tail(colnames(liver2));tail(colnames(heart2)) # Sample IDs are in the right order

# Change the first column name to "Reaction ID" to match the website
# http://marea4galaxy.cloud.ba.infn.it/galaxy
#colnames(heart2)[1] <- "Reaction ID"
#colnames(liver2)[1] <- "Reaction ID"
#head(colnames(liver2));head(colnames(heart2)) # Sample IDs are in the right order

# Round is failing; values appear to not be encoded as floats
lapply(heart2, class) # every column is a character vector

# Convert numeric values from character to float 
dim(heart2[,2:ncol(heart2)]) # 4741 148
class(heart2) # data.table data.frame
heart3 <- sapply(heart2[,2:ncol(heart2)], as.numeric)
lapply(heart3, class)[1:2] # every column is a character vector
heart3 <- as.data.frame(heart3)
class(heart3)
heart3[1:5,1:5]

dim(liver2[,2:ncol(liver2)]) # 4741 148
class(liver2) # data.table data.frame
liver3 <- sapply(liver2[,2:ncol(liver2)], as.numeric)
lapply(liver3, class)[1:2] # every column is a character vector
liver3 <- as.data.frame(liver3)
class(liver3)
liver3[1:5,1:5]

# Append the reaction IDs back onto the dataframe
# Spaces in variables is bad form, but I am doing this to match the example on
# the MaREA for galaxy server because it does not appear to recognize the reactions column
# based on the error message its generating
heart3$"Reaction ID" <- heart2$Reactions
liver3$"Reaction ID" <- liver2$Reactions
# Move "Reactions" column to the front
heart3 <- heart3 %>% select("Reaction ID", everything())
heart3[1:5,1:5]
liver3 <- liver3 %>% select("Reaction ID", everything())
liver3[1:5,1:5]
dim(heart3)==dim(liver3) # TRUE TRUE

# Write to file
#write.table(liver3, "/scratch/mjpete11/mitochondrial_paralogs/linear_models/data/data/liver_RAS_recon22_no_header_complete_cases_only.tsv", quote=FALSE, sep='\t', col.names = NA)
#write.table(heart3, "/scratch/mjpete11/mitochondrial_paralogs/linear_models/data/data/heart_RAS_recon22_no_header_complete_cases_only.tsv", quote=FALSE, sep='\t', col.names = NA)

# Write to file
write.table(liver3, "/scratch/mjpete11/mitochondrial_paralogs/linear_models/data/data/liver_RAS_recon22_complete_cases_only.tsv", quote=FALSE, sep='\t', col.names = TRUE, row.names = FALSE)
write.table(heart3, "/scratch/mjpete11/mitochondrial_paralogs/linear_models/data/data/heart_RAS_recon22_complete_cases_only.tsv", quote=FALSE, sep='\t', col.names = TRUE, row.names = FALSE) 

# Write a small subset to file
small_heart <- heart3[1:10,1:10]
small_liver <- liver3[1:10,1:10]
small_heart[1:5,1:5] 
small_liver[1:5,1:5]

# Write to file
write.table(small_liver, "/scratch/mjpete11/mitochondrial_paralogs/linear_models/data/data/small_liver_RAS_recon22_complete_cases_only.tsv", quote=FALSE, sep='\t', col.names = TRUE, row.names = FALSE) 
write.table(small_heart, "/scratch/mjpete11/mitochondrial_paralogs/linear_models/data/data/small_heart_RAS_recon22_complete_cases_only.tsv", quote=FALSE, sep='\t', col.names = TRUE, row.names = FALSE)
