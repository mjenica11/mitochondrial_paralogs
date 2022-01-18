library(polyester)
library(Biostrings)

# Read in FASTA files of SLC25 mRNA
fasta <- readDNAStringSet("/scratch/mjpete11/simulations/combined.fasta")
head(fasta)

# Fold change matrix; no differential expression between conditions
fold_changes <- matrix(1, 53, 2)
head(fold_changes)

# Reads per transcript
# Control for bias to mapping towards longer reads by including
# the transcript length so the coverage is the same
# ~20x coverage --> reads per transcript = transcriptLength/readLength * 20
readspertx <- round(20 * width(fasta)/100)

# Write number of reads to file
write.csv(readspertx, "readspertx.csv")

# Simulation call
simulate_experiment("/scratch/mjpete11/simulations/combined.fasta",
					readlen=100, paired=TRUE, distr='normal', fraglen=250,
					fragsd=25, error_rate=0.005, bias='none',
					strand_specific=FALSE, set.seed=142,
					num_reps=c(500,500),reads_per_transcript=readspertx, 
					fold_changes=fold_changes, 
					outdir='1000_simulations')

