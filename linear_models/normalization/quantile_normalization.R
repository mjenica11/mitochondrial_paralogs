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

# Read in GTEx counts
counts <- fread("/scratch/mjpete11/linear_models/data/combat_batch_adjusted_counts_zero_filtered.csv", sep=",")

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
setdiff(SLC, counts$'Description') 

#_______________________________________________________________________________ 
# Perform quantile normalization 
#_______________________________________________________________________________ 
rownames(counts) <- counts$V1

quantile_normalisation <- function(df){
		 # Determine the ranks of each column from lowest to highest
		 df_rank <- apply(df,2,rank,ties.method="min")
	     # Sort the original matrix from lowest to highest
 		 df_sorted <- data.frame(apply(df, 2, sort))
		 # Calculate the means
	     df_mean <- apply(df_sorted, 1, mean)
		 # Function to substitute the means into the ranked matrix
         index_to_mean <- function(my_index, my_mean){
			    return(my_mean[my_index])
	  }
      # Substitute means into ranked matrix
	  df_final <- apply(df_rank, 2, index_to_mean, my_mean=df_mean)
	  rownames(df_final) <- rownames(df)
	  return(df_final)
}

new_data <- quantile_normalisation(counts[,3:ncol(counts)])
new_data[1:5,1:5]

# Add gene names column back
new_data <- as.data.frame(new_data)
new_data$"Description" <- counts$"Description"

# Move the gene names column to the front
new_data <- new_data %>% select("Description", everything())

# Write to file
write.csv(new_data, "/scratch/mjpete11/linear_models/data/quantile_normalized_counts_zero_filtered.csv")
