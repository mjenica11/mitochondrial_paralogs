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
manifest <- read.csv("/scratch/mjpete11/mitochondrial_paralogs/linear_models/data/data/sample.tsv", header=TRUE, sep = "\t")

# Make dataframe with sample id, tissue type
# All of the rin number and ischemic time values were missing...
df1 <- data.frame(manifest$"dbgap_sample_id", manifest$"tissue_type")

# Remove all rows with NA in either columns
df2 <- df1[complete.cases(df1), ]

# Read in GTEx counts
counts <- fread("/scratch/mjpete11/mitochondrial_paralogs/linear_models/data/data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct", sep="\t")
counts[1:5,1:5]

# Rename the columns to be more descrptive
colnames(counts)[1] <- "Ensembl_ID"
colnames(counts)[2] <- "Hugo_ID"
counts[1:5,1:5]

# Make a small subset of the unadjusted counts and write to file to tedt with MaREA
smoll <- counts[1:20,1:20]a
# Drop the ensembl_ID column
smoll$Ensembl_ID <- NULL
smoll[1:5,1:5]
# Write small subset to file
write.csv(smoll, "/scratch/mjpete11/mitochondrial_paralogs/linear_models/data/data/small_untransformed_GTEx_counts.csv")

#_______________________________________________________________________________ 
# Perform quantile normalization 
#_______________________________________________________________________________ 
rownames(counts) <- counts$Name
counts[1:5,1:5]

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

# Append the Ensembl ID and Hugo ID name columns back on
new_data$ensembl_ID <- counts$Name
new_data$hugo_ID <- counts$Description

# Make a small subset of the quantile normalized counts and write to file to tedt with MaREA
small <- new_data[1:20,1:20]
write.csv(small, "/scratch/mjpete11/mitochondrial_paralogs/linear_models/data/data/small_quantile_normalized_counts.csv")

# Write to file
write.csv(new_data, "/scratch/mjpete11/mitochondrial_paralogs/linear_models/data/data/quantile_normalized_counts.csv")
