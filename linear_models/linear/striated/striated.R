# Is the proportion of zeros relatively equal across heart and liver samples?

# load libraries
library(tidyr)
library(dplyr)
library(data.table)
library(stringr)

# Read in the total filtered counts matrix 
counts <- fread("/scratch/mjpete11/linear_models/data/filtered_counts.csv", sep=",")

# Read in the organs dataframe because it has the heart and liver sample and subject IDs 
organs <- fread("/scratch/mjpete11/linear_models/data/organs.csv", sep=",") # float

# Subset counts that are from heart or liver samples only
samples <- organs$SAMPID
class(samples) # character vector
heart_liver_counts <- counts[,..samples]
ncol(heart_liver_counts)==length(samples) # TRUE

# Drop duplicate columns
# Identify duplicate columns
dup_cols <- colnames(heart_liver_counts)[duplicated(heart_liver_counts)]

# Read in sample attributes files
metadata <- read.csv("/scratch/mjpete11/linear_models/data/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt", sep="\t")

# Subset heart left ventricle and liver samples into two separate dfs
heart <- metadata$SAMPID[metadata$SMTSD %in% "Heart - Left Ventricle"]
liver <- metadata$SAMPID[metadata$SMTSD %in% "Liver"]

# Which heart sample IDs are present in the colnames of the heart_liver_counts df?
# Instantiate a vector with only those present so the dataframe can be subset
heart_IDs <- intersect(heart, colnames(heart_liver_counts))
liver_IDs <- intersect(liver, colnames(heart_liver_counts))

# Split heart and liver counts dataframe into separate dfs by organ
# Then add the organ variable as a column
heart_df <- subset(heart_liver_counts, select=heart_IDs)
dim(heart_df) # 56200 148
liver_df <- subset(heart_liver_counts, select=liver_IDs)
dim(liver_df) # 56200 148

# Order the columns in the heart and liver dfs so that the paired samples
# have the same index
heart_samples <- heart_IDs[order(match(heart_IDs, liver_IDs))]
head(heart_samples)
head(liver_IDs)
tail(heart_samples)
tail(liver_IDs)

liver_samples <- liver_IDs[order(match(liver_IDs, heart_IDs))]
head(liver_samples)
head(heart_IDs)
tail(liver_samples)
tail(heart_IDs)

# Extract the subject ID from the sample ID 
heart_subjID <- str_sub(heart_samples, start=1L, end=10L)
liver_subjID <- str_sub(liver_samples, start=1L, end=10L)

# Check if they are in the same order
identical(heart_subjID, liver_subjID)

# Re-order the dataframe columns based on the order of the 'samples' vector
setDT(heart_df)
heart_df0 <- setcolorder(heart_df, heart_samples)
head(colnames(heart_df0))
head(heart_samples)

setDT(liver_df)
liver_df0 <- setcolorder(liver_df, liver_samples)
head(colnames(liver_df0))
head(liver_samples)

# Transpose the dataframe
heart_df1 <- t(heart_df0)
liver_df1 <- t(liver_df0)

# Check that the dimentions are the same
dim(heart_df1)==dim(liver_df1) # TRUE TRUE

# Create a copy of each matrix, but replace counts with "success" if value
# is > 0 and "fail" if value is < 0
# Do for one gene and try apply prop.test()
heart_gene <- ifelse(heart_df1[,56200] > 0, "success", "fail")
liver_gene <- ifelse(liver_df1[,56200] > 0, "success", "fail")

# Sometimes you will only get "success" and the table() step will fail
# so I will artificially change the last value to "fail"
# to make sure the function continues
#heart_gene[[56200]] <- "fail"
#liver_gene[[56200]] <- "fail"

# Convert vectors to factors
heart_gene <- as.factor(heart_gene)
liver_gene <- as.factor(liver_gene)

# Make a table of failure and success counts for each gene
tab <- table(heart_gene, liver_gene)
tab

# Proportion test
#test <- prop.test(tab)
#test$p.value # 1 --> proportion of zeros for gene 1 are not 
# sig diff betwen paired heart and liver samples

# Switch to Fisher's exact test because chi-squared test (used by prop.test()
# will fail if there are < 5 instances in either success/fail category
test <- fisher.test(tab)
test$p.value

# Function just to extract the tables
table_function <- function(INDEX){
	heart_gene_x <- ifelse(heart_df1[,INDEX] > 0, "success", "fail")
	liver_gene_x <- ifelse(liver_df1[,INDEX] > 0, "success", "fail")

	# Sometimes you will only get "success" and the table() step will fail
	# because it requires a 2x2 matrix
	# so I will artificially change the last value to "fail" and second to last
	# to "success" to make sure the function continues
	heart_gene_x[[56200]] <- "fail"
	liver_gene_x[[56200]] <- "fail"
	heart_gene_x[[56199]] <- "success"
	liver_gene_x[[56199]] <- "success"

	# Convert vectors to factors
	heart_gene_x <- as.factor(heart_gene_x)
	liver_gene_x <- as.factor(liver_gene_x)

	# Make a table of failure and success counts for each gene
	tab_x <- table(heart_gene_x, liver_gene_x)
	head(tab_x)
	return(tab_x)
}
table_lst <- Map(table_function, INDEX=1:56200)
table_lst

rm(table_lst)
rm(table_function)

# Check if the values are all identical
length(unique(table_lst))==1 # FALSE

# Check if the dimensions are the same for each item
all(sapply(table_lst, function(x) dim(x)==2))==TRUE # TRUE

# Repeat for every gene in the paired heart and liver samples
chi_squared_function <- function(INDEX){
	heart_gene_x <- ifelse(heart_df1[,INDEX] > 0, "success", "fail")
	liver_gene_x <- ifelse(liver_df1[,INDEX] > 0, "success", "fail")

	# Sometimes you will only get "success" and the table() step will fail
	# so I will artificially change the last value to "fail"
	# to make sure the function continues
	heart_gene_x[[56200]] <- "fail"
	liver_gene_x[[56200]] <- "fail"
	heart_gene_x[[56199]] <- "success"
	liver_gene_x[[56199]] <- "success"

	# Convert vectors to factors
	heart_gene_x <- as.factor(heart_gene_x)
	liver_gene_x <- as.factor(liver_gene_x)

	# Make a table of failure and success counts for each gene
	tab_x <- table(heart_gene_x, liver_gene_x)
	head(tab_x)

	# Proportion test
	result <- prop.test(tab_x)
	return(result$p.value)  
}
p_value_lst <- Map(chi_squared_function, INDEX=seq(1, ncol(heart_df), by=1))
p_value_lst

# Check if the values are all identical
length(unique(p_value_lst))==1 # FALSE

rm(chi_squared_function)
rm(p_value_lst)

# Adjust for multiple testing
tests <- length(p_value_lst)

adjusted_p_vals <- lapply(p_value_lst, function(x) x/tests)

# Convert list of p-values to data.frame 
p_val_df <- data.frame(matrix(unlist(adjusted_p_vals), nrow=4, byrow=TRUE),stringsAsFactors=FALSE)
dim(p_val_df)
p_val_df  

# Write to file
write.csv(p_val_df, "/scratch/mjpete11/linear_models/data/two_var_chi_squared_zero_inflated.csv")

