# Combine the simulated RNA-seq expression batches into one dataset

# Load libraries
library(data.table)
library(dplyr)

# Read in simulated batch data
batch1 <- fread("/scratch/mjpete11/mitochondrial_paralogs/linear_models/data/data/simulated_data_salmon_count_matrix_1000_samples.tsv")
batch2 <- fread("/scratch/mjpete11/mitochondrial_paralogs/linear_models/data/data/simulated_batch2_error2_count_matrix_1000_samples.tsv")

head(batch1)
head(batch2)
tail(colnames(batch1))
head(colnames(batch1))
head(colnames(batch2))
tail(colnames(batch2))

head(batch1$V1) # first column has the list of ensembl_IDs
head(batch2$V1) # same

dim(batch1) # 61386 1001
dim(batch2) # 62269 1001

# I have no idea why there are different number of mapped genes
# There should be less genes mapped in batch2::error_rate=2
# because the error rate was super high in the third simulated dataset
# I am guessing there are more mismapped reads and hence more genes falsely detected???
# Drop the decimal points in the ensembl_IDs
batch1$V1 <- sub("\\.\\d+$", "", batch1$V1)
head(batch1$V1) 
batch2$V1 <- sub("\\.\\d+$", "", batch2$V1)
head(batch2$V1)
length(which(batch1$V1 %in% batch2$V1)) # 61358 ensembl gene IDs overlap

# Keep only the rows that are in batch1 in the batch2 dataset
#`%!in%` = Negate(`%in%`)
batch2 <- batch2[(batch2$V1 %in% batch1$V1),]

dim(batch2) # 61358 1001
dim(batch1) # 61386 1001

# Keep only the rows that are in batch2 in the batch1 dataset
# Fingers crossed this does not filter out any SLC25 genes
batch1 <- batch1[(batch1$V1 %in% batch2$V1),]

dim(batch2) # 61358 1001
dim(batch1) # 61358 1001

# Check if the ensembl_IDs equal each other in each dataframe
all(batch1$V1==batch2$V1) # TRUE

# Rename ensembl_ID name column
batch1 <- batch1 %>% rename("ensembl_ID"="V1")
batch2 <- batch2 %>% rename("ensembl_ID"="V1")

head(batch1$ensembl_ID)
head(batch2$ensembl_ID)

# Make sure the dataframe dimensions are equal
dim(batch1)==dim(batch2) # TRUE TRUE

# Order ensembl_IDs to be the same order in both dataframes for merging
batch1 <- batch1 %>% arrange(ensembl_ID, batch2$ensembl_ID)

# Combines into one dataset
combined <- merge(batch1, batch2, by="ensembl_ID", all=TRUE, suffixes=c("_batch1", "_batch2"))
head(combined)
any(is.na(combined)) # FALSE
head(combined$ensembl_ID)
head(colnames(combined))
tail(colnames(combined))
dim(combined) # 61358 2001
combined[1:5,1:5]

# Write to file
write.csv(combined, "/scratch/mjpete11/mitochondrial_paralogs/linear_models/data/data/combined_simulated_batch1_batch2_error2.csv")
