# Batch correction with peer

# Load libraries
library(data.table)
library("peer",lib.loc="/home/mjpete11/miniconda3/envs/peer/lib/R/library/")
#library(peer)

expr <- fread("/scratch/mjpete11/linear_models/data/filtered_counts.csv", sep=",", drop=2)

expr[1:5,1:5]
