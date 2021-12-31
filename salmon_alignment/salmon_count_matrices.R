# Create count matrix of gene and transcript estimates of simulated SLC25 family

# Output
GENES <- "/scratch/mjpete11/salmon_alignment/gene_count_matrix.tsv"
TRANSCRIPTS <- "/scratch/mjpete11/salmon_alignment/transcript_count_matrix.tsv"

# Load packages
library(tximport)
library(GenomicFeatures)

# Vector of sample names
SAMPLES <- sprintf("sample_%02d", 01:20)

# Set path to quant.sf files and check that all files are there
files <- file.path("/scratch/mjpete11/salmon_alignment/mapped_reads", SAMPLES,
				   "quant.sf")
head(files)
all(file.exists(files))

# Make tx2gene table to map gene IDs to transcript IDs
TxDb <- makeTxDbFromGFF(file = "/scratch/mjpete11/salmon_alignment/gencode.v39.annotation.gtf")
k <- keys(TxDb, keytype = "TXNAME")
tx2gene <- select(TxDb, k, "GENEID", "TXNAME") 
head(tx2gene)

# Import salmon output files
# txOut = FALSE to get the list of genes and txOut = TRUE to get the isoforms
txi_gene <- tximport(files, type = "salmon", tx2gene = tx2gene, txOut = FALSE)
txi_trans <- tximport(files, type = "salmon", tx2gene = tx2gene, txOut = TRUE)

counts_gene <- txi_gene$counts
df_gene <- data.frame(counts_gene, check.names = FALSE)

counts_trans <- txi_trans$counts
df_trans <- data.frame(counts_trans, check.names = FALSE)

# Write to tsv
write.table(df_gene, file = GENES, sep = "\t", row.names = TRUE, quote = FALSE)
write.table(df_trans, file = TRANSCRIPTS, sep = "\t", row.names = TRUE, quote = FALSE)
