# Batch correction with peer

# Load libraries
library(data.table)
library("peer",lib.loc="/home/mjpete11/miniconda3/envs/peer/lib/R/library/")
#library(peer)

expr <- fread("/scratch/mjpete11/linear_models/data/filtered_counts.csv", sep=",", drop=2)
expr[1:5,1:5]

dim(expr)

# Create the model object
model=PEER()

# Set the observed data
# NULL response means no error
PEER_setPhenoMean(model, as.matrix(expr))

# In fer k=15 hidden confounders since this is what GTEx guidelines recommend
# for <= 150 samples and there are 148 paired heart and liver samples
PEER_setNk(model, 15)

# Print the number of hidden confounders
PEER_getNk(model)

# Perform the inference
PEER_update(model)

# Figure out what the model object returned is so I can write it to file
print(class(model))
