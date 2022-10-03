# Perform quantile normalization on count data

# Libraries
library(data.table)
library(tidyverse)
library(reshape)
library(sva)

#_______________________________________________________________________________ 
# Make matrix of genes, samples, counts, organs, and batch effects
#_______________________________________________________________________________ 
# Read in GTEx manifest
manifest <- read.csv("/scratch/mjpete11/linear_models/data/sample.tsv", header=TRUE, sep = "\t")

# Make dataframe with sample id, tissue type
# All of the rin number and ischemic time values were missing...
df1 <- data.frame(manifest$"dbgap_sample_id", manifest$"tissue_type")

# Remove all rows with NA in either columns
df2 <- df1[complete.cases(df1), ]

# Read in GTEx counts
counts <- fread("/scratch/mjpete11/linear_models/data/combat_batch_adjusted_counts.csv", sep=",")

#_______________________________________________________________________________ 
# Perform quantile normalization 
#_______________________________________________________________________________ 
rownames(counts) <- counts$Name

quantile_normalisation <- function(df){
		 df_rank <- apply(df,2,rank,ties.method="min")
 		 df_sorted <- data.frame(apply(df, 2, sort))
	     df_mean <- apply(df_sorted, 1, mean)

         index_to_mean <- function(my_index, my_mean){
			    return(my_mean[my_index])
	  }

	  df_final <- apply(df_rank, 2, index_to_mean, my_mean=df_mean)
	  rownames(df_final) <- rownames(df)
	  return(df_final)
}

#new_data <- quantile_normalisation(counts[,2:ncol(counts)])
#new_data <- counts[1:5,3:7]
new_data <- quantile_normalisation(counts[,3:ncol(counts)])
new_data[1:5,1:5]

# Write to file
write.csv(new_data, "/scratch/mjpete11/linear_models/data/quantile_normalized_counts.csv")
