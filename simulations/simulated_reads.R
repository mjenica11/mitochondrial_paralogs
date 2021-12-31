library(polyester)
library(Biostrings)

# Read in FASTA files of SLC25 mRNA
fasta <- readDNAStringSet("/scratch/mjpete11/simulations/combined.fasta")
head(fasta)

# Fold change matrix; no differential expression between conditions
fold_changes <- matrix(1, 53, 2)
head(fold_changes)

# Simulation call
simulate_experiment("/scratch/mjpete11/simulations/combined.fasta",
					readlen=100, paired=TRUE, distr='normal', fraglen=250,
					fragsd=25, error <- rate=0.005, bias='none',
					strand_specific=FALSE, meanmodel=TRUE, num_reps=c(10,10),
					fold_changes=fold_changes, outdir='simulated_reads1')
