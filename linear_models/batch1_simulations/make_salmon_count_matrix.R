# Create count matrix from salmon results of the simulated data. 
# Note: if complaint about select(), unload dplyr
setwd("/scratch/mjpete11/linear_models/batch1_simulations/")

# Load packages                                                                 
library(tximport)                                                               
library(GenomicFeatures) 
library(AnnotationDbi)
library(edgeR)   
library(readr)

# Input
METADATA = "/scratch/mjpete11/linear_models/data/filtered_counts_organ_metadata.csv"

# Output
GENE_FILE = "/scratch/mjpete11/linear_models/data/simulated_data_salmon_count_matrix.tsv"

# Read Metadata CSV.                                                            
Samples = read.csv(METADATA, header = TRUE)

# Set rownames of metadata object equal to sample names.                        
rownames(Samples) <- Samples$Sample   

# Set path to quant.sf files and check that all files are there.                
genes <- c("SLC25A1", "SLC25A2", "SLC25A3", "SLC25A4", "SLC25A5", "SLC25A6", "SLC25A7", "SLC25A8", "SLC25A9", "SLC25A10", "SLC25A11", "SLC25A12", "SLC25A13", "SLC25A14", "SLC25A15", "SLC25A16", "SLC25A17", "SLC25A18", "SLC25A19", "SLC25A20", "SLC25A21", "SLC25A22", "SLC25A23", "SLC25A24", "SLC25A25", "SLC25A26", "SLC25A27", "SLC25A28", "SLC25A29", "SLC25A30", "SLC25A31", "SLC25A32", "SLC25A33", "SLC25A34", "SLC25A35", "SLC25A36", "SLC25A37", "SLC25A38", "SLC25A39", "SLC25A40", "SLC25A41", "SLC25A42", "SLC25A43", "SLC25A44", "SLC25A45", "SLC25A46", "SLC25A47", "SLC25A48", "SLC25A49", "SLC25A50", "SLC25A51", "SLC25A52", "SLC25A53")
files <- file.path("/scratch/mjpete11/linear_models/batch1_simulations/individual_simulations", genes, "quant.sf")
head(files)                                                                     
all(file.exists(files)) # TRUE

# Set the names of the file paths object equal to the sample names.             
names(files) <- Samples$genes                                                  

# Make tx2gene table to map gene IDs to transcript IDs.                         
TxDb <- makeTxDbFromGFF(file = "/scratch/mjpete11/linear_models/batch1_simulations/gencode.v41.annotation.gff3.gz")
k <- keys(TxDb, keytype = "TXNAME")                                             
tx2gene <- select(TxDb, k, "GENEID", "TXNAME") # write output with gene ID/ trans ID (if one exists)                                 
head(tx2gene)                                                                   

# Import salmon output files.
# If you want the list of genes, add: txOut = FALSE. To get list of isoforms, add: txOut = TRUE
Txi_Gene <- tximport(files, type = "salmon", tx2gene = tx2gene, txOut = FALSE, ignoreAfterBar=TRUE)      
Counts_Gene <- Txi_Gene$counts
df_Gene <- data.frame(Counts_Gene, check.names = FALSE)

# Write to TSV
write.table(df_Gene, file = GENE_FILE, sep = "\t", row.names = TRUE, quote = FALSE)
