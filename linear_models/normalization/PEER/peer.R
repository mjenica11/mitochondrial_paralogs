# Batch correction with peer on only the 53 SLC25 genes in the paired heart and liver samples
# Tutorial: https://github.com/PMBio/peer/wiki/Tutorial

# Load libraries
library(data.table)
library("peer",lib.loc="/home/mjpete11/miniconda3/envs/peer/lib/R/library/")

# Read in filtered counts
#expr <- fread("/scratch/mjpete11/linear_models/data/filtered_counts.csv", sep=",", drop=2)
#expr[1:5,1:5]

# Drop the index column
#expr$V1 <- NULL
#dim(expr)

# Write a small portion to file for Liu
#small <- expr[1:20,1:20]
#write.csv(small, "/scratch/mjpete11/linear_models/data/sample_batch_voom_qnorm_matrix.csv")

counts <- fread("/scratch/mjpete11/linear_models/data/filtered_counts.csv", sep=",")

# Drop the index column
counts$V1 <- NULL

# Generate the design matrix
# Samples are rows and columns are covariates of interest (heart and liver)

# List of SLC25 paralogs
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

# Subset the SLC25 genes
counts <- as.data.frame(counts)
sub_df <- counts[counts$"Description" %in% SLC, ]
dim(sub_df)

# Read in count df of just SLC genes
#sub_df <- read.csv("/scratch/mjpete11/linear <- models/data/SLC <- df <- voom <- combat <- seq.csv", sep=",")
names(sub_df) <- gsub(x=names(sub_df), pattern="\\.", replacement="-")

# All present!
setdiff(SLC, sub_df$'Description') 

# Use this file to create the design matrix
organs <- read.csv("/scratch/mjpete11/linear_models/data/filtered_counts_organ_metadata.csv", sep=",")

# Drop samples from the counts df that are not from paired heart or liver
keeps <- as.vector(organs$SAMPID)
expr <- subset(sub_df, select=keeps) 

# Are there the same number of samples in the counts df and the organs metadata df? 
ncol(expr)==nrow(organs) # TRUE
dim(expr)
dim(organs)

# Subset to just the covariates and convert to a matrix
covs <- as.matrix(organs[, c("SMRIN", "SMTSISCH")])

# Create the model object
model=PEER()

# Set the observed data
# NULL response means no error
PEER_setPhenoMean(model, as.matrix(expr))
print("setting the phenotype mean")

# Include known covariates
# This sets the first C factors to be fixed to the obseved covariates and extends
# the hideen factor matrix X to have additional C columns
#PEER_setCovariates(model, as.matrix(covs))

# In fer k=15 hidden confounders since this is what GTEx guidelines recommend
# for <= 150 samples and there are 148 paired heart and liver samples
PEER_setNk(model, 30)
print("setNK")

# Print the number of hidden confounders
PEER_getNk(model)
print("getNK")

# Perform the inference
PEER_update(model)
print("update model")

# Figure out what the model object returned is so I can write it to file
print("printing class of model object")
class(model) # p_PEER_cSPARSEFA (??? what on earth is this object...)
#write.csv(model, "/scratch/mjpete11/linear_models/data/peer_model.csv")

# Get the posterior mean of the inferred confounders (N x K matrix)
factors <- PEER_getX(model)
print("printing class of factors object")
class(factors) # matrix
write.csv(factors, "/scratch/mjpete11/linear_models/data/peer_factors.csv")
print("Wrote factors")

# Get the posterior weights of the inferred confounders (G x K) matrix
weights_obj <- PEER_getW(model)
print("printing class of weights object")
class(weights_obj) # matrix
write.csv(weights_obj, "/scratch/mjpete11/linear_models/data/peer_weights.csv")
print("wrote weights object")

# Get the precision (inverse variance) of the weights (K x 1) matrix
precision <- PEER_getAlpha(model)
print("printing class of precision object")
class(precision) # matrix
write.csv(precision, "/scratch/mjpete11/linear_models/data/peer_precision.csv")
print("wrote class of precision object")

# Get the residual dataset (N x G) matrix
residuals <- PEER_getResiduals(model)
print("printing residuals")
class(residuals) # matrix
write.csv(residuals, "/scratch/mjpete11/linear_models/data/peer_residuals.csv")
print("wrote residuals to file")

# Plot 
PEER_plotModel(model)
