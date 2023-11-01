library(polyester)
library(Biostrings)

# Read in FASTA files of SLC25 mRNA
fasta <- readDNAStringSet("/scratch/mjpete11/linear_models/batch2_simulations/combined.fasta")
head(fasta)

# Fold change matrix; no differential expression between conditions
fold_changes <- matrix(1, 53, 2) 
head(fold_changes)
dim(fold_changes)

# Reads per transcript
# Control for bias to mapping towards longer reads by including
# the transcript length so the coverage is the same
# ~20x coverage --> reads per transcript = transcriptLength/readLength * 20
readspertx <- round(20 * width(fasta)/100)
length(readspertx)
tmp <- matrix(c(nrow=readspertx, ncol=readspertx))
dim(tmp)
dim(tmp)==dim(fold_changes)

# Write number of reads to file
write.csv(readspertx, "readspertx.csv")

# Simulation call
simulate_experiment("/scratch/mjpete11/linear_models/batch2_simulations/combined.fasta",
					readlen=100, paired=TRUE, distr='normal', fraglen=250,
					fragsd=25, error_rate=0.05, add_platform_error='illumina4', bias='none',
					strand_specific=FALSE, set.seed=97, # change seed compared to batch 1
					num_reps=c(500,500),reads_per_transcript=readspertx, 
					fold_changes=fold_changes, 
					outdir='1000_samples')

