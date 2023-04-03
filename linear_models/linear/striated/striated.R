# Are the striated samples always the same samples?

# load libraries
library(tidyr)
library(dplyr)
library(data.table)

# Are the striated genes correlating with those samples?
zero_inflated_genes <- c("SLC25A2", "UCP1", "SLC25A21", "SLC2531", 
						   "SLC25A41", "SLC25A47", "SLC25A48", "SLC25A52")
not_zero_inflated_genes <- c("SLC25A1", "SLC25A3", "SLC25A4", "SLC25A5", 
							   "SLC25A6", "UCP2", "SLC25UCP3", "SLC25A10", 
							   "SLC25A11", "SLC25A12", "SLC25A13", "SLC25A14", 
							   "SLC25A15", "SLC25A16", "SLC25A17", "SLC25A18", 
							   "SLC25A19", "SLC25A20", "SLC25A22", "SLC25A23", 
							   "SLC25A24", "SLC25A25", "SLC25A26", "SLC25A27", 
							   "SLC25A28", "SLC25A29", "SLC25A30", "SLC25A32",
				  			   "SLC25A33", "SLC25A34", "SLC25A35", "SLC25A36",
				  			   "SLC25A37", "SLC25A38", "SLC25A39", "SLC25A40",
				  			   "SLC25A42", "SLC25A43", "SLC25A44", "SLC25A45", 
							   "SLC25A46", "SLC25A49", "SLC25A50", "SLC25A51", "SLC25A53")

# Read in voom quantile normalized counts
organs <- fread("/scratch/mjpete11/linear_models/data/voom_qnorm_counts1.csv", sep=",") # float

# List of all samples
all_samples <- colnames(filtered_counts)[3:ncol(filtered_counts)]

# Values of zero-inflated SLCs 
zero_inflated_counts <- subset(organs, subset=gene  %in% zero_inflated_genes)
zero_inflated_gene_correlations <- cor(x=zero_inflated_counts$counts, y=zero_inflated_counts$SAMPID)

not_zero_inflated_gene_correlations <- cor(x=not_zero_inflated_genes, y=all_samples)

#________________________________________________________________________________ 
# I realized my sparsity matrix idea doesn't make sense because it's specific 
# genes that are zero-inflated, not the samples so it doesn't make sense
# to compare the samples...keeping the following code in case it ends up being useful
#________________________________________________________________________________ 

# Read in striated samples after blocking by batch with voom and
# striated samples without batch correction
# Write distribution of gene values to file
# These are the genes from the non-batch corrected heart and liver samples
# that had an excess of zero values (i.e. the "striated" samples)
no_batch_A52 <- read.csv("/scratch/mjpete11/linear_models/data/sparsity_matrix_files/A52_no_batch.csv", sep=",")
no_batch_UCP1 <- read.csv("/scratch/mjpete11/linear_models/data/sparsity_matrix_files/UCP1_no_batch.csv", sep=",")
no_batch_A31 <- read.csv("/scratch/mjpete11/linear_models/data/sparsity_matrix_files/A31_no_batch.csv", sep=",")
no_batch_A47 <- read.csv("/scratch/mjpete11/linear_models/data/sparsity_matrix_files/A47_no_batch.csv", sep=",")
no_batch_A2 <- read.csv("/scratch/mjpete11/linear_models/data/sparsity_matrix_files/A2_no_batch.csv", sep=",")
no_batch_A48 <- read.csv("/scratch/mjpete11/linear_models/data/sparsity_matrix_files/A48_no_batch.csv", sep=",")
no_batch_A21 <- read.csv("/scratch/mjpete11/linear_models/data/sparsity_matrix_files/A21_no_batch.csv", sep=",")
no_batch_A41 <- read.csv("/scratch/mjpete11/linear_models/data/sparsity_matrix_files/A41_no_batch.csv", sep=",")

# These are the genes from the non-batch corrected heart and liver samples
# that did not have an excess of zero values (i.e. the "non-striated" samples)
no_batch_A1 <- read.csv("/scratch/mjpete11/linear_models/data/sparsity_matrix_files/A1_no_batch.csv", sep=",")
no_batch_A3 <- read.csv("/scratch/mjpete11/linear_models/data/sparsity_matrix_files/A3_no_batch.csv", sep=",")
no_batch_A4 <- read.csv("/scratch/mjpete11/linear_models/data/sparsity_matrix_files/A4_no_batch.csv", sep=",")
no_batch_A5 <- read.csv("/scratch/mjpete11/linear_models/data/sparsity_matrix_files/A5_no_batch.csv", sep=",")
no_batch_A6 <- read.csv("/scratch/mjpete11/linear_models/data/sparsity_matrix_files/A6_no_batch.csv", sep=",")
no_batch_UCP2 <- read.csv("/scratch/mjpete11/linear_models/data/sparsity_matrix_files/UCP2_no_batch.csv", sep=",")
no_batch_UCP3 <- read.csv("/scratch/mjpete11/linear_models/data/sparsity_matrix_files/UCP3_no_batch.csv", sep=",")
no_batch_A10 <- read.csv("/scratch/mjpete11/linear_models/data/sparsity_matrix_files/A10_no_batch.csv", sep=",")
no_batch_A11 <- read.csv("/scratch/mjpete11/linear_models/data/sparsity_matrix_files/A11_no_batch.csv", sep=",")
no_batch_A12 <- read.csv("/scratch/mjpete11/linear_models/data/sparsity_matrix_files/A12_no_batch.csv", sep=",")
no_batch_A13 <- read.csv("/scratch/mjpete11/linear_models/data/sparsity_matrix_files/A13_no_batch.csv", sep=",")
no_batch_A14 <- read.csv("/scratch/mjpete11/linear_models/data/sparsity_matrix_files/A14_no_batch.csv", sep=",")
no_batch_A15 <- read.csv("/scratch/mjpete11/linear_models/data/sparsity_matrix_files/A15_no_batch.csv", sep=",")
no_batch_A16 <- read.csv("/scratch/mjpete11/linear_models/data/sparsity_matrix_files/A16_no_batch.csv", sep=",")
no_batch_A17 <- read.csv("/scratch/mjpete11/linear_models/data/sparsity_matrix_files/A17_no_batch.csv", sep=",")
no_batch_A18 <- read.csv("/scratch/mjpete11/linear_models/data/sparsity_matrix_files/A18_no_batch.csv", sep=",")
no_batch_A19 <- read.csv("/scratch/mjpete11/linear_models/data/sparsity_matrix_files/A19_no_batch.csv", sep=",")
no_batch_A20 <- read.csv("/scratch/mjpete11/linear_models/data/sparsity_matrix_files/A20_no_batch.csv", sep=",")
no_batch_A22 <- read.csv("/scratch/mjpete11/linear_models/data/sparsity_matrix_files/A22_no_batch.csv", sep=",")
no_batch_A23 <- read.csv("/scratch/mjpete11/linear_models/data/sparsity_matrix_files/A23_no_batch.csv", sep=",")
no_batch_A24 <- read.csv("/scratch/mjpete11/linear_models/data/sparsity_matrix_files/A24_no_batch.csv", sep=",")
no_batch_A25 <- read.csv("/scratch/mjpete11/linear_models/data/sparsity_matrix_files/A25_no_batch.csv", sep=",")
no_batch_A26 <- read.csv("/scratch/mjpete11/linear_models/data/sparsity_matrix_files/A26_no_batch.csv", sep=",")
no_batch_A27 <- read.csv("/scratch/mjpete11/linear_models/data/sparsity_matrix_files/A27_no_batch.csv", sep=",")
no_batch_A28 <- read.csv("/scratch/mjpete11/linear_models/data/sparsity_matrix_files/A28_no_batch.csv", sep=",")
no_batch_A29 <- read.csv("/scratch/mjpete11/linear_models/data/sparsity_matrix_files/A29_no_batch.csv", sep=",")
no_batch_A29 <- read.csv("/scratch/mjpete11/linear_models/data/sparsity_matrix_files/A29_no_batch.csv", sep=",")
no_batch_A30 <- read.csv("/scratch/mjpete11/linear_models/data/sparsity_matrix_files/A30_no_batch.csv", sep=",")
no_batch_A31 <- read.csv("/scratch/mjpete11/linear_models/data/sparsity_matrix_files/A31_no_batch.csv", sep=",")
no_batch_A32 <- read.csv("/scratch/mjpete11/linear_models/data/sparsity_matrix_files/A32_no_batch.csv", sep=",")
no_batch_A33 <- read.csv("/scratch/mjpete11/linear_models/data/sparsity_matrix_files/A33_no_batch.csv", sep=",")
no_batch_A34 <- read.csv("/scratch/mjpete11/linear_models/data/sparsity_matrix_files/A34_no_batch.csv", sep=",")
no_batch_A35 <- read.csv("/scratch/mjpete11/linear_models/data/sparsity_matrix_files/A35_no_batch.csv", sep=",")
no_batch_A36 <- read.csv("/scratch/mjpete11/linear_models/data/sparsity_matrix_files/A36_no_batch.csv", sep=",")
no_batch_A37 <- read.csv("/scratch/mjpete11/linear_models/data/sparsity_matrix_files/A37_no_batch.csv", sep=",")
no_batch_A38 <- read.csv("/scratch/mjpete11/linear_models/data/sparsity_matrix_files/A38_no_batch.csv", sep=",")
no_batch_A39 <- read.csv("/scratch/mjpete11/linear_models/data/sparsity_matrix_files/A39_no_batch.csv", sep=",")
no_batch_A40 <- read.csv("/scratch/mjpete11/linear_models/data/sparsity_matrix_files/A40_no_batch.csv", sep=",")
no_batch_A42 <- read.csv("/scratch/mjpete11/linear_models/data/sparsity_matrix_files/A42_no_batch.csv", sep=",")
no_batch_A43 <- read.csv("/scratch/mjpete11/linear_models/data/sparsity_matrix_files/A43_no_batch.csv", sep=",")
no_batch_A44 <- read.csv("/scratch/mjpete11/linear_models/data/sparsity_matrix_files/A44_no_batch.csv", sep=",")
no_batch_A45 <- read.csv("/scratch/mjpete11/linear_models/data/sparsity_matrix_files/A45_no_batch.csv", sep=",")
no_batch_A46 <- read.csv("/scratch/mjpete11/linear_models/data/sparsity_matrix_files/A46_no_batch.csv", sep=",")
no_batch_A49 <- read.csv("/scratch/mjpete11/linear_models/data/sparsity_matrix_files/A49_no_batch.csv", sep=",")
no_batch_A50 <- read.csv("/scratch/mjpete11/linear_models/data/sparsity_matrix_files/A50_no_batch.csv", sep=",")
no_batch_A51 <- read.csv("/scratch/mjpete11/linear_models/data/sparsity_matrix_files/A51_no_batch.csv", sep=",")
no_batch_A53 <- read.csv("/scratch/mjpete11/linear_models/data/sparsity_matrix_files/A53_no_batch.csv", sep=",")

# Write distribution of gene values to file
# These are the genes from the voom adjusted and blocked by batch heart and liver samples
# that had an excess of zero values (i.e. the "striated" samples)
batch_UCP1 <- read.csv("/scratch/mjpete11/linear_models/data/sparsity_matrix_files/UCP1_batch.csv", sep=",")
batch_A52 <- read.csv("/scratch/mjpete11/linear_models/data/sparsity_matrix_files/A52_batch.csv", sep=",")
batch_A31 <- read.csv("/scratch/mjpete11/linear_models/data/sparsity_matrix_files/A31_batch.csv", sep=",")
batch_A2 <- read.csv("/scratch/mjpete11/linear_models/data/sparsity_matrix_files/A2_batch.csv", sep=",")
batch_A47 <- read.csv("/scratch/mjpete11/linear_models/data/sparsity_matrix_files/A47_batch.csv", sep=",")
batch_A48 <- read.csv("/scratch/mjpete11/linear_models/data/sparsity_matrix_files/A48_batch.csv", sep=",")
batch_A21 <- read.csv("/scratch/mjpete11/linear_models/data/sparsity_matrix_files/A21_batch.csv", sep=",")
batch_A41 <- read.csv("/scratch/mjpete11/linear_models/data/sparsity_matrix_files/A41_batch.csv", sep=",")

# These are the genes from the non-batch corrected heart and liver samples
# that did not have an excess of zero values (i.e. the "non-striated" samples)
batch_A1 <- read.csv("/scratch/mjpete11/linear_models/data/sparsity_matrix_files/A1_batch.csv", sep=",")
batch_A3 <- read.csv("/scratch/mjpete11/linear_models/data/sparsity_matrix_files/A3_batch.csv", sep=",")
batch_A4 <- read.csv("/scratch/mjpete11/linear_models/data/sparsity_matrix_files/A4_batch.csv", sep=",")
batch_A5 <- read.csv("/scratch/mjpete11/linear_models/data/sparsity_matrix_files/A5_batch.csv", sep=",")
batch_A6 <- read.csv("/scratch/mjpete11/linear_models/data/sparsity_matrix_files/A6_batch.csv", sep=",")
batch_UCP2 <- read.csv("/scratch/mjpete11/linear_models/data/sparsity_matrix_files/UCP2_batch.csv", sep=",")
batch_UCP3 <- read.csv("/scratch/mjpete11/linear_models/data/sparsity_matrix_files/UCP3_batch.csv", sep=",")
batch_A10 <- read.csv("/scratch/mjpete11/linear_models/data/sparsity_matrix_files/A10_batch.csv", sep=",")
batch_A11 <- read.csv("/scratch/mjpete11/linear_models/data/sparsity_matrix_files/A11_batch.csv", sep=",")
batch_A12 <- read.csv("/scratch/mjpete11/linear_models/data/sparsity_matrix_files/A12_batch.csv", sep=",")
batch_A13 <- read.csv("/scratch/mjpete11/linear_models/data/sparsity_matrix_files/A13_batch.csv", sep=",")
batch_A14 <- read.csv("/scratch/mjpete11/linear_models/data/sparsity_matrix_files/A14_batch.csv", sep=",")
batch_A15 <- read.csv("/scratch/mjpete11/linear_models/data/sparsity_matrix_files/A15_batch.csv", sep=",")
batch_A16 <- read.csv("/scratch/mjpete11/linear_models/data/sparsity_matrix_files/A16_batch.csv", sep=",")
batch_A17 <- read.csv("/scratch/mjpete11/linear_models/data/sparsity_matrix_files/A17_batch.csv", sep=",")
batch_A18 <- read.csv("/scratch/mjpete11/linear_models/data/sparsity_matrix_files/A18_batch.csv", sep=",")
batch_A19 <- read.csv("/scratch/mjpete11/linear_models/data/sparsity_matrix_files/A19_batch.csv", sep=",")
batch_A20 <- read.csv("/scratch/mjpete11/linear_models/data/sparsity_matrix_files/A20_batch.csv", sep=",")
batch_A22 <- read.csv("/scratch/mjpete11/linear_models/data/sparsity_matrix_files/A22_batch.csv", sep=",")
batch_A23 <- read.csv("/scratch/mjpete11/linear_models/data/sparsity_matrix_files/A23_batch.csv", sep=",")
batch_A24 <- read.csv("/scratch/mjpete11/linear_models/data/sparsity_matrix_files/A24_batch.csv", sep=",")
batch_A25 <- read.csv("/scratch/mjpete11/linear_models/data/sparsity_matrix_files/A25_batch.csv", sep=",")
batch_A26 <- read.csv("/scratch/mjpete11/linear_models/data/sparsity_matrix_files/A26_batch.csv", sep=",")
batch_A27 <- read.csv("/scratch/mjpete11/linear_models/data/sparsity_matrix_files/A27_batch.csv", sep=",")
batch_A28 <- read.csv("/scratch/mjpete11/linear_models/data/sparsity_matrix_files/A28_batch.csv", sep=",")
batch_A29 <- read.csv("/scratch/mjpete11/linear_models/data/sparsity_matrix_files/A29_batch.csv", sep=",")
batch_A29 <- read.csv("/scratch/mjpete11/linear_models/data/sparsity_matrix_files/A29_batch.csv", sep=",")
batch_A30 <- read.csv("/scratch/mjpete11/linear_models/data/sparsity_matrix_files/A30_batch.csv", sep=",")
batch_A31 <- read.csv("/scratch/mjpete11/linear_models/data/sparsity_matrix_files/A31_batch.csv", sep=",")
batch_A32 <- read.csv("/scratch/mjpete11/linear_models/data/sparsity_matrix_files/A32_batch.csv", sep=",")
batch_A33 <- read.csv("/scratch/mjpete11/linear_models/data/sparsity_matrix_files/A33_batch.csv", sep=",")
batch_A34 <- read.csv("/scratch/mjpete11/linear_models/data/sparsity_matrix_files/A34_batch.csv", sep=",")
batch_A35 <- read.csv("/scratch/mjpete11/linear_models/data/sparsity_matrix_files/A35_batch.csv", sep=",")
batch_A36 <- read.csv("/scratch/mjpete11/linear_models/data/sparsity_matrix_files/A36_batch.csv", sep=",")
batch_A37 <- read.csv("/scratch/mjpete11/linear_models/data/sparsity_matrix_files/A37_batch.csv", sep=",")
batch_A38 <- read.csv("/scratch/mjpete11/linear_models/data/sparsity_matrix_files/A38_batch.csv", sep=",")
batch_A39 <- read.csv("/scratch/mjpete11/linear_models/data/sparsity_matrix_files/A39_batch.csv", sep=",")
batch_A40 <- read.csv("/scratch/mjpete11/linear_models/data/sparsity_matrix_files/A40_batch.csv", sep=",")
batch_A42 <- read.csv("/scratch/mjpete11/linear_models/data/sparsity_matrix_files/A42_batch.csv", sep=",")
batch_A43 <- read.csv("/scratch/mjpete11/linear_models/data/sparsity_matrix_files/A43_batch.csv", sep=",")
batch_A44 <- read.csv("/scratch/mjpete11/linear_models/data/sparsity_matrix_files/A44_batch.csv", sep=",")
batch_A45 <- read.csv("/scratch/mjpete11/linear_models/data/sparsity_matrix_files/A45_batch.csv", sep=",")
batch_A46 <- read.csv("/scratch/mjpete11/linear_models/data/sparsity_matrix_files/A46_batch.csv", sep=",")
batch_A49 <- read.csv("/scratch/mjpete11/linear_models/data/sparsity_matrix_files/A49_batch.csv", sep=",")
batch_A50 <- read.csv("/scratch/mjpete11/linear_models/data/sparsity_matrix_files/A50_batch.csv", sep=",")
batch_A51 <- read.csv("/scratch/mjpete11/linear_models/data/sparsity_matrix_files/A51_batch.csv", sep=",")
batch_A53 <- read.csv("/scratch/mjpete11/linear_models/data/sparsity_matrix_files/A53_batch.csv", sep=",")



# List object to contain dataframes
# list of dataframes of samples that were not batch corrected
no_batch_lst <- list(no_batch_A1, no_batch_A2, no_batch_A3, no_batch_A4,
					 no_batch_A5, no_batch_A6, no_batch_UCP1, no_batch_UCP2,
					 no_batch_UCP3, no_batch_A10, no_batch_A11, no_batch_A12,
					 no_batch_A13, no_batch_A14, no_batch_A15, no_batch_A16,
					 no_batch_A17, no_batch_A18, no_batch_A19, no_batch_A20,
					 no_batch_A21, no_batch_A22, no_batch_A23, no_batch_A24,
					 no_batch_A25, no_batch_A26, no_batch_A27, no_batch_A28,
					 no_batch_A29, no_batch_A30, no_batch_A31, no_batch_A32,
					 no_batch_A33, no_batch_A34, no_batch_A35, no_batch_A36,
					 no_batch_A37, no_batch_A38, no_batch_A39, no_batch_A40,
					 no_batch_A41, no_batch_A42, no_batch_A43, no_batch_A44,
					 no_batch_A45, no_batch_A46, no_batch_A47, no_batch_A48,
					 no_batch_A49, no_batch_A50, no_batch_A51, no_batch_A52, no_batch_A53)

# List of dataframes of samples that were corrected by blocking by batch with voom
batch_lst <- list(batch_A1, batch_A2, batch_A3, batch_A4,
				  batch_A5, batch_A6, batch_UCP1, batch_UCP2,
				  batch_UCP3, batch_A10, batch_A11, batch_A12,
				  batch_A13, batch_A14, batch_A15, batch_A16,
				  batch_A17, batch_A18, batch_A19, batch_A20,
				  batch_A21, batch_A22, batch_A23, batch_A24,
				  batch_A25, batch_A26, batch_A27, batch_A28,
				  batch_A29, batch_A30, batch_A31, batch_A32,
				  batch_A33, batch_A34, batch_A35, batch_A36,
				  batch_A37, batch_A38, batch_A39, batch_A40,
				  batch_A41, batch_A42, batch_A43, batch_A44,
				  batch_A45, batch_A46, batch_A47, batch_A48,
				  batch_A49, batch_A50, batch_A51, batch_A52, batch_A53)

# Function to drop the 'value' column in the samples that were not corrected for
# batch effects so they can be merged with the batch corrected dataframes
drop_col <- function(dat){
		dat$value <- NULL
		return(dat)
}
no_batch_lst1 <- Map(drop_col, no_batch_lst)

# Combine gene distribution dataframes into a list object
# Striated (zero-inflated) samples: A52, UCP1, A31, A47, A2, A48, A21, A41
# These are present in both the batch and non-batch corrected samples
striated_lst <- list(no_batch_lst1[[52]], no_batch_lst1[[7]], no_batch_lst1[[31]], 
					 no_batch_lst1[[47]], no_batch_lst1[[2]], no_batch_lst1[[48]],
					 no_batch_lst1[[21]], no_batch_lst1[[41]], batch_lst[[52]], 
					 batch_lst[[7]], batch_lst[[31]], batch_lst[[2]], 
					 batch_lst[[47]], batch_lst[[48]], batch_lst[[21]], batch_lst[[41]]) 

# Combine list of zero inflated dfs into one dataframe
striated_res <- Reduce(function(x, y) merge(x, y, all = TRUE), striated_lst)

# Find the intersection of the subject ID calls
striated_res0 <- striated_res[duplicated(striated_res$SAMPID),]
head(striated_res0)

# Are there any samples that are present in one condition (e.g. no batch correction)
# vs another (e.g. blocking for batch by voom)
#res2 <- !res0$SAMPID %in% duplicated(res0$SAMPID) 
#lst1 <- list(no_batch_A52, no_batch_UCP1, no_batch_A31, no_batch_A47, 
#			no_batch_A2, no_batch_A48, no_batch_A21, no_batch_A41)

#lst2 <- list(batch_A52, batch_UCP1, batch_A31, batch_A2,
#			batch_A47, batch_A48, batch_A21, batch_A41) 

#res_no_batch <- Reduce(function(x, y) merge(x, y, all = TRUE), lst1)
#res_batch <- Reduce(function(x, y) merge(x, y, all = TRUE), lst2)

#samples <- unique(res_batch$SAMPID)

#any(!res_no_batch$SAMPID %in% samples)==TRUE # FALSE

# Sort the df rows by the sample IDs
#res1 <- res0[order(res0$SAMPID),]
#head(res1)
#tail(res1)

# Write to file
#write.csv(res1, "/scratch/mjpete11/linear_models/data/striated_intersection.csv", sep=",")

# Read in file
samples <- read.csv("/scratch/mjpete11/linear_models/data/striated_intersection.csv", sep=",")

# Read in filtered counts (genes with >5 counts in at least one sample)
counts <- fread("/scratch/mjpete11/linear_models/data/filtered_counts.csv", sep=",")

# Calculate the sparsity score for the non-striated (not zero-inflated) samples
# sparsity = count zero elements / total elements
# Test
striated_samples <- unique(striated_res0$SAMPID)
s1 <- striated_samples[1]
sample1 <- counts[[s1]] 
sparse1 <- sum(counts[[s1]]==0)/length(counts[[s1]])

# Scale up
sparsity <- function(SAMPLE_ID) {
		res <- sum(counts[[SAMPLE_ID]]==0)/length(counts[[SAMPLE_ID]]) 
		return(res)
}
sparsity_lst <- Map(sparsity, SAMPLE_ID=striated_samples)

# Range
range(sparsity_lst) # 0.43 to 0.60

# Add an extra element so I have a multiple and can convert from a list to a matrix object
sparsity_lst0 <- append(sparsity_lst, "NA", after=length(sparsity_lst))

# Reshape into a matrix
mat <- matrix(unlist(sparsity_lst0), ncol=8, byrow=TRUE)

# Write to file
write.csv(mat, "/scratch/mjpete11/linear_models/data/sparsity_matrix.csv")

# Make a list of the non-striated (not zero-inflated) samples 
not_striated_lst <- list(no_batch_lst1[[1]], no_batch_lst1[[3]], no_batch_lst1[[4]], 
					     no_batch_lst1[[5]], no_batch_lst1[[6]], no_batch_lst1[[8]], 
           				 no_batch_lst1[[9]], no_batch_lst1[[10]], no_batch_lst1[[11]], 
						 no_batch_lst1[[12]], no_batch_lst1[[13]], no_batch_lst1[[14]], 
						 no_batch_lst1[[15]], no_batch_lst1[[16]], no_batch_lst1[[17]], 
						 no_batch_lst1[[18]], no_batch_lst1[[19]], no_batch_lst1[[20]], 
						 no_batch_lst1[[22]], no_batch_lst1[[23]], no_batch_lst1[[24]], 
						 no_batch_lst1[[25]], no_batch_lst1[[26]], no_batch_lst1[[27]], 
						 no_batch_lst1[[28]], no_batch_lst1[[29]], no_batch_lst1[[30]], 
						 no_batch_lst1[[32]], no_batch_lst1[[33]], no_batch_lst1[[34]], 
						 no_batch_lst1[[35]], no_batch_lst1[[36]], no_batch_lst1[[37]], 
						 no_batch_lst1[[38]], no_batch_lst1[[39]], no_batch_lst1[[40]],  
						 no_batch_lst1[[42]], no_batch_lst1[[43]], no_batch_lst1[[44]], 
						 no_batch_lst1[[45]], no_batch_lst1[[46]], no_batch_lst1[[49]], 
						 no_batch_lst1[[50]], no_batch_lst1[[51]], no_batch_lst1[[53]], 
						 batch_lst[[1]], batch_lst[[3]], batch_lst[[4]], 
					     batch_lst[[5]], batch_lst[[6]], batch_lst[[8]], 
           				 batch_lst[[9]], batch_lst[[10]], batch_lst[[11]], 
						 batch_lst[[12]], batch_lst[[13]], batch_lst[[14]], 
						 batch_lst[[15]], batch_lst[[16]], batch_lst[[17]], 
						 batch_lst[[18]], batch_lst[[19]], batch_lst[[20]], 
						 batch_lst[[22]], batch_lst[[23]], batch_lst[[24]], 
						 batch_lst[[25]], batch_lst[[26]], batch_lst[[27]], 
						 batch_lst[[28]], batch_lst[[29]], batch_lst[[30]], 
						 batch_lst[[32]], batch_lst[[33]], batch_lst[[34]], 
						 batch_lst[[35]], batch_lst[[36]], batch_lst[[37]], 
						 batch_lst[[38]], batch_lst[[39]], batch_lst[[40]],  
						 batch_lst[[42]], batch_lst[[43]], batch_lst[[44]], 
						 batch_lst[[45]], batch_lst[[46]], batch_lst[[49]], 
						 batch_lst[[50]], batch_lst[[51]],  batch_lst[[53]])

# Combine list of dfs into one dataframe
not_striated_res <- Reduce(function(x, y) merge(x, y, all = TRUE), not_striated_lst)

# Make the sparsity matrix of the non-striated samples (samples with extra zeros)
# so I can compare them with the non-zero inflated samples
non_striated_samples <- unique(not_striated_res$SAMPID) 

# Create a sparsity matrix out of the non-striated samples
non_sparsity_lst <- Map(sparsity, SAMPLE_ID=non_striated_samples)

# Range
range(non_sparsity_lst) # 0.43 to 0.60; same as the non-zero inflated samples...

# Reshape into a matrix
non_sparse_mat <- matrix(unlist(non_sparsity_lst), ncol=8, byrow=TRUE)

# Write to file
write.csv(non_sparse_mat, "/scratch/mjpete11/linear_models/data/non_striated_samples_sparsity_matrix.csv")

