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
#METADATA = "/scratch/mjpete11/linear_models/data/filtered_counts_organ_metadata.csv"

# Output
GENE_FILE = "/scratch/mjpete11/linear_models/data/simulated_data_salmon_count_matrix_1000_samples.tsv"

# Read Metadata CSV.                                                            
#Samples = read.csv(METADATA, header = TRUE)

# Create a vector of sample names
numbers <- seq(1, 1000)
for (i in 1:10){
		if (i <= 9) {
				numbers[i] <- paste0("0",i)
		} else {
				numbers[i] <- i
		}
}
strings <- paste("sample", numbers, sep="_")
head(strings);tail(strings)

# Set rownames of metadata object equal to sample names.                        
#rownames(Samples) <- Samples$Sample   

# Set path to quant.sf files and check that all files are there.                
files <- file.path("/scratch/mjpete11/linear_models/batch1_simulations/individual_simulations", strings, "quant.sf")
head(files);tail(files)
all(file.exists(files)) # TRUE

# Set the names of the file paths object equal to the sample names.             
names(files) <- strings                                                  

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
