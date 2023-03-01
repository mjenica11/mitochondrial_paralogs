# Are the striated samples always the same samples?

# Read in striated samples after blocking by batch with voom and
# striated samples without batch correction
# Write distribution of gene values to file
no_batch_A52 <- read.csv("/scratch/mjpete11/linear_models/data/striated_A52_no_batch.csv", sep=",")
no_batch_UCP1 <- read.csv("/scratch/mjpete11/linear_models/data/striated_UCP1_no_batch.csv", sep=",")
no_batch_A31 <- read.csv("/scratch/mjpete11/linear_models/data/striated_A31_no_batch.csv", sep=",")
no_batch_A47 <- read.csv("/scratch/mjpete11/linear_models/data/striated_A47_no_batch.csv", sep=",")
no_batch_A2 <- read.csv("/scratch/mjpete11/linear_models/data/striated_A2_no_batch.csv", sep=",")
no_batch_A48 <- read.csv("/scratch/mjpete11/linear_models/data/striated_A48_no_batch.csv", sep=",")
no_batch_A21 <- read.csv("/scratch/mjpete11/linear_models/data/striated_A21_no_batch.csv", sep=",")
no_batch_A41 <- read.csv("/scratch/mjpete11/linear_models/data/striated_A41_no_batch.csv", sep=",")

# Write distribution of gene values to file
batch_UCP1 <- read.csv("/scratch/mjpete11/linear_models/data/striated_UCP1_batch_voom.csv", sep=",")
batch_A52 <- read.csv("/scratch/mjpete11/linear_models/data/striated_A52_batch_voom.csv", sep=",")
batch_A31 <- read.csv("/scratch/mjpete11/linear_models/data/striated_A31_batch_voom.csv", sep=",")
batch_A2 <- read.csv("/scratch/mjpete11/linear_models/data/striated_A2_batch_voom.csv", sep=",")
batch_A47 <- read.csv("/scratch/mjpete11/linear_models/data/striated_A47_batch_voom.csv", sep=",")
batch_A48 <- read.csv("/scratch/mjpete11/linear_models/data/striated_A48_batch_voom.csv", sep=",")
batch_A21 <- read.csv("/scratch/mjpete11/linear_models/data/striated_A21_batch_voom.csv", sep=",")
batch_A41 <- read.csv("/scratch/mjpete11/linear_models/data/striated_A41_batch_voom.csv", sep=",")

# Drop the column that the no_batch dfs have but the batch dfs don't have
no_batch_UCP1$value <- NULL
no_batch_A52$value <- NULL
no_batch_A31$value <- NULL
no_batch_A2$value <- NULL
no_batch_A47$value <- NULL
no_batch_A48$value <- NULL
no_batch_A21$value <- NULL
no_batch_A41$value <- NULL

# Combine gene distribution dataframes into a list object
lst <- list(no_batch_A52, no_batch_UCP1, no_batch_A31, no_batch_A47, 
			no_batch_A2, no_batch_A48, no_batch_A21, no_batch_A41,
			batch_A52, batch_UCP1, batch_A31, batch_A2,
			batch_A47, batch_A48, batch_A21, batch_A41) 

# Combine list of dfs into one dataframe
res <- Reduce(function(x, y) merge(x, y, all = TRUE), lst)

# Find the intersection of the subject ID calls
res0 <- res[duplicated(res$SAMPID),]
head(res0)

# Sort the df rows by the sample IDs
rm(res1)
res1 <- res0[order(res0$SAMPID),]
head(res1)
tail(res1)

# Write to file
write.csv(res1, "/scratch/mjpete11/linear_models/data/striated_intersection.csv", sep=",")
