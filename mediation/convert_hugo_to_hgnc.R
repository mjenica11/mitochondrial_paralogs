# Convert the HUGO gene names column to the HUGO gene ID column for use with MaREA

# Load the libraries
library(data.table)
library(hgnc)
library(dplyr)

# Read in the small heart and liver count data
heart <- fread("/scratch/mjpete11/mitochondrial_paralogs/linear_models/data/data/GTEx_heart_filtered_counts.tsv", sep="\t")
liver <- fread("/scratch/mjpete11/mitochondrial_paralogs/linear_models/data/data/GTEx_liver_filtered_counts.tsv", sep="\t")

# Read in the gene names
hugo_name_vec <- heart$Hugo_ID
head(hugo_name_vec);tail(hugo_name_vec) 
length(hugo_name_vec) # 51259
class(hugo_name_vec) # list 

# Import the HGNC dataset
hgnc_data <- import_hgnc_dataset()
class(hgnc_data) # tbl_df tbl data.frame
class(hgnc_data$hgnc_id) #character 
head(hgnc_data)
names(hgnc_data)
head(hgnc_data$hgnc_id) # This is what MaREA is designed to accept as input
head(hgnc_data$name)
head(hgnc_data$symbol)
length(hgnc_data$name) # 43842 ... there are less genes in the HGNC dataset than in my GTEx dataset ...

# Check that the vector of HUGO IDs in my GTEx normalized count data is present in the HGNC database
any(hugo_name_vec %in% hgnc_data$name) == TRUE # TRUE 
any(hugo_name_vec %in% hgnc_data$symbol) == TRUE # TRUE 
all(hugo_name_vec %in% hgnc_data$name) == TRUE # FALSE 
all(hugo_name_vec %in% hgnc_data$symbol) == TRUE # FALSE 

# How many genes are present in the GTEx unfiltered counts and in
# in the hgnc database if we search using the HUGO gene name
length(intersect(hugo_name_vec, hgnc_data$symbol)) # 31395

# Check if the MCF family are in the HGNC database as I have the names recorded
SLC <- c("SLC25A1", "SLC25A2", "SLC25A3", "SLC25A4", "SLC25A5", "SLC25A6", 
"UCP1", "UCP2", "UCP3", "SLC25A10", "SLC25A11", "SLC25A12", 
"SLC25A13", "SLC25A14", "SLC25A15", "SLC25A16", "SLC25A17", 
"SLC25A18", "SLC25A19", "SLC25A20", "SLC25A21", "SLC25A22", 
"SLC25A23", "SLC25A24", "SLC25A25", "SLC25A26", "SLC25A27", 
"SLC25A28", "SLC25A29", "SLC25A30", "SLC25A31", "SLC25A32", 
"SLC25A33", "SLC25A34", "SLC25A35", "SLC25A36", "SLC25A37", 
"SLC25A38", "SLC25A39", "SLC25A40", "SLC25A41", "SLC25A42", 
"SLC25A43", "SLC25A44", "SLC25A45", "SLC25A46", "SLC25A47",
"SLC25A48", "MTCH1", "MTCH2", "SLC25A51", "SLC25A52", "SLC25A53")

all(SLC %in% hgnc_data$symbol) == TRUE # FALSE 
which(SLC %in% hgnc_data$symbol) == TRUE # FALSE 

# Function to check if MCF genes are in the HGNC database
check_MCF <- function(x){
    boolean <- SLC[x] %in% hgnc_data$symbol
    boolean
}
check_MCF(x=7)
mapply(check_MCF, seq(1:53)) # The MCF family is in the HUGO database as I've been recording them 

# Create a logical vector indicating if the HUGO genes in the GTEx dataframe is present in the HGNC database 
hgnc_id_subset <- hugo_name_vec %in% hgnc_data$symbol
any(hgnc_id_subset)==TRUE # TRUE
hgnc_id_subset

# Subset using logical vectors
# This results in the HGNC names that correspond to the HUGO IDs
hgnc_name_subset <- hugo_name_vec[hgnc_id_subset]
head(hgnc_name_subset);tail(hgnc_name_subset)
class(hgnc_name_subset) # character

# Subset the HGNC IDs that correspond to the HGNC names
subset_hgnc_df <- hgnc_data[hgnc_data$symbol %in% hgnc_name_subset,]
dim(subset_hgnc_df) # 31395 55 ; why are there so few matches???
subset_hgnc_df[1:5,1:5]

# Drop the rows from the GTEx dataframe that are not present in the HUGO database
# I have to drop them because there will be no HUGO ID and then I can't use it with MaREA
# First, rename the HUGO gene name column to 'symbol' so I can inner_join() with a matching column name
heart_subset_df<- heart[heart$Hugo_ID %in% subset_hgnc_df$symbol,]
heart_subset_df[1:5,1:5]
dim(heart_subset_df) # 31568 431

liver_subset_df<- liver[liver$Hugo_ID %in% subset_hgnc_df$symbol,]
liver_subset_df[1:5,1:5]
dim(liver_subset_df) # 31568 277

# Rename the HUGO gene name column to 'name' so I can bind the Hugo_ID column from the subset_hgnc_df
colnames(liver_subset_df)[1] <- "symbol"
colnames(heart_subset_df)[1] <- "symbol"

# Append a column with the actual Hugo IDs (I originally named the column incorrectly)
head(subset_hgnc_df$symbol)
head(subset_hgnc_df$hgnc_id)
hugoID_subset_heart <- left_join(heart_subset_df, subset_hgnc_df[,c("hgnc_id", "symbol")], by="symbol")
ncol(hugoID_subset_heart)==ncol(heart_subset_df)+1 # TRUE
# Move the hugoID column to the front
tail(colnames(hugoID_subset_heart))
hugoID_subset_heart <- hugoID_subset_heart %>% select(hgnc_id, everything())
hugoID_subset_heart[1:5,1:5]

# Check that the order of the HGNC names match the order of the HGNC name vector
# The order needs to match or else you will append the wrong HGNC name/ID to the wrong entry!
subset_hgnc_df$name==hgnc_name_subset # Mostly FALSE

# Reorder the hgnc dataframe rows to be in the same order as the hgnc_name_subset vector
reordered_subset_df <- subset_hgnc_df %>% slice(order(factor(name, levels=hgnc_name_subset)))
head(reordered_subset_df$name)
all(reordered_subset_df$name==hgnc_name_subset) # TRUE

# Append the reordered HGNC dataframe to the heart and liver count dataframes
reordered_subset_df[1:5,]
reordered_subset_df$name
reordered_subset_df$hgnc_id
hgnc_name_subset
heart$HGNC_ID <- reordered_subset_df$hgnc_id
liver$HGNC_ID <- reordered_subset_df$hgnc_id
heart[1:5,1:5]
liver[1:5,1:5]

# Move the HGNC ID column to the front of each dataframe
heart %>% relocate(HGNC_ID)
