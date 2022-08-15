# Perform combat normalization on count data

# Libraries
library(data.table)
library(tidyverse)
library(reshape)
library(sva)

# Read in GTEx counts
counts <- fread("/scratch/mjpete11/linear_models/data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct", sep="\t")

# Read in metadata with batch variables
manifest <- read.csv("/scratch/mjpete11/linear_models/data/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt", header=TRUE, sep = "\t")

# Drop samples from manifest that are missing in the count data
rownames(manifest) <- manifest$SAMPID
manifest <- manifest[(manifest$SAMPID %in% colnames(counts)[3:ncol(counts)]),]
samples <- manifest$SAMPID

# Only keep samples that are present in the manifest 
counts <- counts[,..samples]

# Make a list of known batch variables; genotype or expression batch ID 
batch <- manifest$SMGEBTCH
#batch <- manifest$SMGEBTCH[1:18]

# Test count df
#test <- counts[,3:20]

# Create a model matrix for the adjustment variables
# Only added and intercept term, since I am only adjusting for batch
modcombat <- model.matrix(~1, data=manifest)
#modcombat <- model.matrix(~1, data=manifest[1:18,])

# Apply ComBat to the data, using parametric empirical Bayesian adjustment
combat_edata <- ComBat(dat=counts, batch=batch, 
					   mod=modcombat, par.prior=TRUE)
#combat_edata <- ComBat(dat=test, batch=batch, 
#					   mod=modcombat, par.prior=TRUE)

# Write to file
write.csv(combat_edata, "/scratch/mjpete11/linear_models/results2/combat_batch_adjusted_counts.csv")
