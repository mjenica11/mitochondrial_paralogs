# Combine the simulated RNA-seq expression batches into one dataset

# Load libraries
library(data.table)
library(dplyr)

# Read in simulated batch data
batch1 <- fread("/scratch/mjpete11/linear_models/data/simulated_data_salmon_count_matrix_1000_samples.tsv")
batch2 <- fread("/scratch/mjpete11/linear_models/data/simulated_batch2_salmon_count_matrix_1000_samples.tsv")

head(batch1)
head(batch2)
tail(colnames(batch1))
head(colnames(batch1))
tail(colnames(batch2))

head(batch1$V1)
head(batch2$gene)

# Rename ensembl_ID name column
batch1 <- batch1 %>% rename("ensembl_ID"="V1")
batch2 <- batch2 %>% rename("ensembl_ID"="gene")

head(batch1$ensembl_ID)
head(batch2$ensembl_ID)

# Make sure the dataframe dimensions are equal
batch2$V1 <- NULL
dim(batch1) # 61386 1001
dim(batch2) # 61386 1001
dim(batch1)==dim(batch2)

# Make sure ensembl_ID rows are in the same order for both dataframes
all(batch1$ensembl_ID==batch2$ensembl_ID) # FALSE

# Order ensembl_IDs to be the same order in both dataframes for merging
batch1 <- batch1 %>% arrange(ensembl_ID, batch2$ensembl_ID)

# Combines into one dataset
combined <- merge(batch1, batch2, by="ensembl_ID", all=TRUE, suffixes=c("_batch1", "_batch2"))
head(combined)
any(is.na(combined)) # FALSE
head(combined$ensembl_ID)
head(colnames(combined))
tail(colnames(combined))
dim(combined) # 61386 2001

# Write to file
write.csv(combined, "/scratch/mjpete11/linear_models/data/combined_simulated_batch1_batch2.csv")
