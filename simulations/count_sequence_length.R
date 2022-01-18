library(Biostrings)

s <- readAAStringSet("/scratch/mjpete11/simulations/combined.fasta", "fasta")
summary(nchar(s))
