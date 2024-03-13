# Convert the HUGO gene names column to the HUGO gene ID column for use with MaREA

# Load the libraries
library(data.table)
library(hgnc)

# Read in the small heart and liver count data
heart <- fread("/scratch/mjpete11/mitochondrial_paralogs/linear_models/data/data/GTEx_heart_filtered_counts.tsv", sep="\t")
liver <- fread("/scratch/mjpete11/mitochondrial_paralogs/linear_models/data/data/GTEx_liver_filtered_counts.tsv", sep="\t")

# Read in the gene names
gene_names <- as.list(heart$Hugo_ID)
gene_names[1:5]
class(gene_names) # list 

# Import the HGNC dataset
hgnc_data <- import_hgnc_dataset()
class(hgnc_data) # tbl_df tbl data.frame
head(hgnc_data)
names(hgnc_data)
head(hgnc_data$hgnc_id) # This is what MaREA is designed to accept as input
head(hgnc_data$name)

# Filter data by gene names
any(gene_names %in% hgnc_data$name) == TRUE # TRUE 
any(hgnc_data$name %in% gene_names) == TRUE
hgnc_filtered <- hgnc_data %>% hgnc_id %in% gene_names
