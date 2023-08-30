# Filter simulated batch expression counts 

# Libraries
library(data.table)
library(tidyverse)
library(reshape)
library(sva)
library(dplyr)

# Read in GTEx counts
counts <- fread("/scratch/mjpete11/linear_models/data/combined_simulated_batch1_batch2.csv", sep=",")
dim(counts) # 61386 2002
counts[1:5,1:5]
counts[1:5,2:5]

# Drop the index column
counts$V1 <- NULL

# Expression filter: Keep a gene if it has >x counts in at least 1 sample 
# Skip the first column because that is just the gene name column
expression_filter <- function(DF, thresh){
				DAT <- DF[rowSums(DF[,2:ncol(DF)]>=thresh, na.rm=TRUE)>0, ]
				return(DAT)
}

# Apply function
counts2 <- expression_filter(DF=counts, thresh=5)
dim(counts2) # 53 2001

# Use another function to test if it is working...
tmp <- counts2[rowSums(counts2[,2:ncol(counts2)]>=5, na.rm=TRUE)>0,] 
dim(tmp) # 53 2001

################################ test code ####################################
# Create a sample dataframe
dat <- data.frame(a = c("X", "Y", "Z", "W"),
				 b = c(1, 2, 3, 4),
				 c = c(5, 4, 3, 2),
				 d = c(9, NA, 11, 12))

# Drop rows with no value greater than or equal to 5 in columns b, c, and d
dat <- dat[rowSums(dat[, 2:ncol(dat)] >= 5, na.rm = TRUE) > 0, ]
dim(dat) #  3 4
# Drop rows with no value greater than or equal to 5
dat1 <- dat[rowSums(dat[, 2:ncol(dat)] >= 5, na.rm = TRUE) > 0, ]
dim(dat1) #  3 4
dat2 <- expression_filter(DF=dat, thresh=5)
dim(dat2) # 3 4
#...Both functions are definitely working...
# Diagnosed the problem! Needed to explicitely skip the gene name column
################################ test code ####################################

# How many genes are left if you set the threshold to 10?
counts_10 <- expression_filter(DF=counts2, thresh=10)
dim(counts_10) # 53 2001
test0 <- counts2[rowSums(counts2[,2:ncol(counts2)]>=10, na.rm=TRUE)>0,] # Still got the same result
dim(test0) # 53 2001.. results can be replicated

# Add hugo_ID column
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

counts2$hugo_ID <- SLC

# Add ensembl_ID without version number column
ensembl_IDs <- c("ENSG00000100075", "ENSG00000120329", "ENSG00000075415",
				"ENSG00000151729", "ENSG00000005022", "ENSG00000169100",
				"ENSG00000109424", "ENSG00000175567", "ENSG00000175564",
				"ENSG00000183048", "ENSG00000108528", "ENSG00000115840",
				"ENSG00000004864", "ENSG00000102078", "ENSG00000102743",
				"ENSG00000122912", "ENSG00000100372", "ENSG00000182902",
				"ENSG00000125454", "ENSG00000178537", "ENSG00000183032",
				"ENSG00000177542", "ENSG00000125648", "ENSG00000085491",
				"ENSG00000148339", "ENSG00000144741", "ENSG00000153291",
				"ENSG00000155287", "ENSG00000197119", "ENSG00000174032",
				"ENSG00000151475", "ENSG00000164933", "ENSG00000171612",
				"ENSG00000162461", "ENSG00000125434", "ENSG00000114120",
				"ENSG00000147454", "ENSG00000144659", "ENSG00000013306",
				"ENSG00000075303", "ENSG00000181240", "ENSG00000181035",
				"ENSG00000077713", "ENSG00000160785", "ENSG00000162241",
				"ENSG00000164209", "ENSG00000140107", "ENSG00000145832",
				"ENSG00000137409", "ENSG00000109919", "ENSG00000122696",
				"ENSG00000141437", "ENSG00000269743")

counts2$ensembl_ID <- ensembl_IDs

# Move the hugo_ID column to the front
counts2 <- counts2 %>% dplyr::select(hugo_ID, everything())
counts2[1:5,1:5]

# Write filtered counts to file
write.csv(counts2, "/scratch/mjpete11/linear_models/data/simulated_filtered_counts.csv")

