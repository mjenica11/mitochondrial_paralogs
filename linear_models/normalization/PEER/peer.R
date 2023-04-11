# Batch correction with peer
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

# Read in voom and quantile normalized and filtered expression data
expr <- fread("/scratch/mjpete11/linear_models/data/batch_voom_qnorm_matrix.csv")
expr$V1 <- NULL
expr[1:5,1:5]
dim(expr)

# Write a small portion to file for Liu
small <- expr[1:20,1:20]
write.csv(small, "/scratch/mjpete11/linear_models/data/sample_batch_voom_qnorm_matrix.csv")

# Read in dataframe containing the covariates
#organs <- fread("/scratch/mjpete11/linear_models/data/voom_qnorm_counts1.csv", sep=",")

# Subset to just the covariates and convert to a matrix
#covs <- as.matrix(organs[, c("SMRIN", "SMTSISCH")])

# Create the model object
model=PEER()

# Set the observed data
# NULL response means no error
PEER_setPhenoMean(model, as.matrix(expr))

# Include known covariates
# This sets the first C factors to be fixed to the obseved covariates and extends
# the hideen factor matrix X to have additional C columns
#PEER_setCovariates(model, as.matrix(covs))

# In fer k=15 hidden confounders since this is what GTEx guidelines recommend
# for <= 150 samples and there are 148 paired heart and liver samples
PEER_setNk(model, 4)

# Print the number of hidden confounders
PEER_getNk(model)

# Perform the inference
PEER_update(model)

# Figure out what the model object returned is so I can write it to file
print("printing class of model object")
class(model) # p_PEER_cSPARSEFA (??? what on earth is this object...)
#write.csv(model, "/scratch/mjpete11/linear_models/data/peer_model.csv")

# Get the posterior mean of the inferred confounders (N x K matrix)
factors <- PEER_getX(model)
print("printing class of factors object")
class(factors) # matrix
write.csv(factors, "/scratch/mjpete11/linear_models/data/peer_factors.csv")

# Get the posterior weights of the inferred confounders (G x K) matrix
weights_obj <- PEER_getW(model)
print("printing class of weights object")
class(weights_obj) # matrix
write.csv(weights_obj, "/scratch/mjpete11/linear_models/data/peer_weights.csv")

# Get the precision (inverse variance) of the weights (K x 1) matrix
precision <- PEER_getAlpha(model)
print("printing class of precision object")
class(precision) # matrix
write.csv(precision, "/scratch/mjpete11/linear_models/data/peer_precision.csv")

# Get the residual dataset (N x G) matrix
residuals <- PEER_getResiduals(model)
print("printing class of precision object")
class(residuals) # matrix
write.csv(residuals, "/scratch/mjpete11/linear_models/data/peer_residuals.csv")

# Plot 
plot(precision)
