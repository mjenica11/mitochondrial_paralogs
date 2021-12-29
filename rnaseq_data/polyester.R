# Simulate expression from ADP/ATP paralogs

#BiocManager::install("polyester")
library(polyester)
library(Biostrings)

# FASTA annotation
fasta_file <- system.file('extdata', 'chr22.fa', package='polyester')
fasta = readDNAStringSet(fasta_file)

# Define fold changes
fold_changes = matrix(c(4,4,rep(1,18),1,1,4,4,rep(1,16)), nrow=20)
head(fold_changes)

# Subset the FASTA file to the first 20 transcripts
small_fasta <- fasta[1:20]
writeXStringSet(small_fasta, 'chr22_small.fa')

# ~10x coverage ---> reads per transcripts = transcriptLength/readLength * 10
# Here all transcripts will have ~equal FPKM
readspertx = round(10 * width(small_fasta) / 100)

# Simulation call
simulate_experiment(fasta='chr22_small.fa', 
					reads_per_transcripts=readspertx, num_reps=c(10,10),
					fold_changes, outdir='simulated_reads')

