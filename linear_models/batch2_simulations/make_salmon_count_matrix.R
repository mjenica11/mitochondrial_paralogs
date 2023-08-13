# Create count matrix from salmon results of the simulated data. 
# Note: if complaint about select(), unload dplyr
setwd("/scratch/mjpete11/linear_models/batch1_simulations/")

# Load packages                                                                 
library(tximport)                                                               
library(GenomicFeatures) 
library(AnnotationDbi)
library(edgeR)   
library(readr)
library(dplyr)
library(data.table)

# Input
#METADATA = "/scratch/mjpete11/linear_models/data/filtered_counts_organ_metadata.csv"

# Output
GENE_FILE = "/scratch/mjpete11/linear_models/data/simulated_batch2_salmon_count_matrix_1000_samples.tsv"
#TEST = "/scratch/mjpete11/linear_models/data/simulated_data_salmon_count_matrix_20_samples.tsv"

# Read Metadata CSV.                                                            
#Samples = read.csv(METADATA, header = TRUE)

# Create a vector of sample names
numbers <- seq(1, 1000)
for (i in 1:1000){
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
files <- file.path("/scratch/mjpete11/linear_models/batch2_simulations/individual_simulations", strings, "quant.sf")
head(files);tail(files)
all(file.exists(files)) # TRUE
which(!file.exists(files)) # TRUE

# Set the names of the file paths object equal to the sample names.             
names(files) <- strings                                                  

# Make tx2gene table to map gene IDs to transcript IDs.                         
TxDb <- makeTxDbFromGFF(file = "/scratch/mjpete11/linear_models/batch2_simulations/gencode.v41.annotation.gff3.gz")
k <- keys(TxDb, keytype = "TXNAME")                                             
tx2gene <- select(TxDb, k, "GENEID", "TXNAME") # write output with gene ID/ trans ID (if one exists)                                 
head(tx2gene)                                                                   

# Import salmon output files.
# If you want the list of genes, add: txOut = FALSE. To get list of isoforms, add: txOut = TRUE
Txi_Gene <- tximport(files, type = "salmon", tx2gene = tx2gene, txOut = FALSE, ignoreAfterBar=TRUE)      
Txi_Gene$counts[1:5,1:5]
Counts_Gene <- Txi_Gene$counts
dim(Counts_Gene)
df_Gene <- data.frame(Counts_Gene, check.names = TRUE)
dim(df_Gene) # 61386 1000

# Name the gene ID column 
#df_Gene$gene <- row.names(df_Gene)
df_Gene1 <- df_Gene
df_Gene1 <- tibble::rownames_to_column(df_Gene1, "ensembl_ID")
df_Gene1[1:5,1:5]
df_Gene1 <- df_Gene1 %>% relocate(gene) %>% head()
dim(df_Gene1)

# Make a smaller test dataframe and save rownames as a column so it can be 
# read back in correctly
rm(test)
df_Gene1$row_name <- row.names(df_Gene1)
test <- df_Gene1[1:20,1:20]
test[1:5,1:5]
dim(df_Gene1) # 61386 1001

object.size("sample_1000") # 120 bytes
object.size("ENSG00000000003.15") # 136 bytes
sizes <- sapply(colnames(df_Gene1), object.size)
median(sizes) # 120 bytes
sizes_y <- sapply(rownames(df_Gene1), object.size)
median(sizes_y) # 136 bytes

# Write to TSV
write.table(df_Gene1, file = GENE_FILE, sep = "\t", row.names = TRUE, quote = FALSE)
#write.table(df_Gene1, file = TEST, sep = "\t", row.names = TRUE, quote = FALSE)

# Try reading the table back in
rm(dat)
dat <- fread("/scratch/mjpete11/linear_models/data/simulated_data_salmon_count_matrix_20_samples.tsv", sep="\t")
dat[1:5,1:5]
