# Make a dataframe of the SLC25 and metabolic biomarkers

# Libraries
library(plyr)
library(limma)
library(data.table)
library(tidyverse)
library(reshape)
library(stringr)
library(ggplot2)
library(ggpubr)
library(edgeR)

# Read in voom + qnorm adjusted filtered GTEx count matrix
voom_matrix <- fread("/scratch/mjpete11/mitochondrial_paralogs/linear_models/data/data/batch_voom_qnorm_matrix.csv", sep=",") 
voom_matrix[1:10,1:10]
dim(voom_matrix) # 51259 298 

# Read in sample attributes files
manifest <- read.csv("/scratch/mjpete11/mitochondrial_paralogs/linear_models/data/data/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt", sep="\t")

# Drop the index column
voom_matrix$V1 <- NULL

# Convert matrix to data.frame
voom_df <- as.data.frame(voom_matrix)
dim(voom_df) # 51259 297

# List of SLC25 paralogs
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

biomarkers <- c("SLC16A1", "SLC16A2", "HK1", "HK2", "HK3", "GPI", "PFKM",
                "PFKP", "ALDOA", "ALDOB", "TPI1", "GAPDH", "PGK1","PGAM1", 
                "PGAM5", "PGAM4", "ENO1","PKLR", "LDHA", "LDHB", "LDHC", 
                "LDHD", "CPT1A","CPT1B", "CPT2", "ACADVL", "ACADL", "ACADM", 
                "ACADS","ACADSB", "ACAD11", "ACAD8", "ACAD10", "ACAD9", 
                "ECHS1","ECH1", "ECI1", "ECI2", "HADHA", "PDHA1", 
                "PDHA2","PDHB", "PDP1", "PDHX", "CS", "ACO1", "ACO2", "IDH3G",
                "IDH1", "IDH2", "IDH3B", "OGDH", "SUCLG2","SUCLG1", "SUCLA2", 
                "SDHB", "SDHD", "SDHA", "SDHC", "SDHAF2","FH", "MDH1", "MDH2", "MT-ND1")

# Double check the vector lengths
length(SLC) # 53
length(biomarkers) # 65

# Add gene name column
gene_names <- fread("/scratch/mjpete11/mitochondrial_paralogs/linear_models/data/data/filtered_counts.csv", sep=",",
					select=c("Description"))
head(gene_names)

# Convert gene_names from df to vector
gene_names <- gene_names[["Description"]] 
class(gene_names)
dim(voom_df) # 51259 297 
length(gene_names) # 51259

# Add the gene names (hugo_IDs) column to the voom adjusted counts dataframe
voom_df$hugo_IDs <- gene_names 

# Combine the metabolic biomarkers and SLC25s into one vector
genes_to_subset <- c(SLC, biomarkers)
genes_to_subset
length(genes_to_subset) # 117

# Subset the SLC25 genes
sub_df <- voom_df[voom_df$hugo_IDs %in% genes_to_subset, ]
sub_df[1:5,1:5]

# Number of genes remaining
dim(voom_df) # 51259 298 
dim(sub_df) # 117 298 

# Which biomarker genes are missing?
setdiff(sub_df$hugo_IDs, genes_to_subset) # 0
which(genes_to_subset %in% unique(sub_df$hugo_IDs)==FALSE)
genes_to_subset[92]
sub_df$hugo_IDs[!sub_df$hugo_IDs %in% gene_names] # 0
length(intersect(sub_df$hugo_IDs, gene_names)) # 117 

# Are there any duplicated row names?
any(duplicated(sub_df$hugo_names))==TRUE # FALSE

# Keep only unique colnames
voom_df <- subset(sub_df, select=unique(colnames(voom_df)))
dim(voom_df) # 117 297 (148 samples * 2 conditions + 1 gene name columns)

# Move the gene names column to the front
voom_df <- voom_df %>% select("hugo_IDs", everything())
voom_df[1:5,1:5]

# Subset heart left ventricle and liver samples into two separate dfs
heart <- manifest[manifest$SMTSD %in% "Heart - Left Ventricle", ]
liver <- manifest[manifest$SMTSD %in% "Liver", ]

# Subset the voom + qnorm adjusted matrix by the sample IDs that match the IDs in file3 df
heart2 <- voom_df %>% select(contains(heart$SAMPID))
liver2 <- voom_df %>% select(contains(liver$SAMPID))

# How many heart/liver samples are there?
ncol(heart2) # 148
ncol(liver2) # 148

dim(heart2) # 117 148
dim(liver2) # 117 148

# Append the gene ensemble ID and common name
heart3 <- cbind('gene'=voom_df$'hugo_IDs', heart2)
liver3 <- cbind('gene'=voom_df$'hugo_IDs', liver2)

# Add a columns with the organ in each df
heart3$organ <- 'heart' 
liver3$organ <- 'liver' 

dim(heart3) # 117 150
dim(liver3) # 117 150

head(heart3)
head(liver3)

# Reshape dataframe so that gene, organ, sample ID and voom adjusted values
# are columns
heart4 <- melt(data = heart3,
			   id.vars = c("gene", "organ"),
			   measure.vars = colnames(heart3)[3:ncol(heart3)-1],
			   variable.name = "samples",
			   value.name = "counts")

liver4 <- melt(data = liver3,
			   id.vars = c("gene", "organ"),
			   measure.vars = colnames(liver3)[3:ncol(liver3)-1],
			   variable.name = "samples",
			   value.name = "counts")

# Change the name of the 'variable' column to 'SAMPID' to match columns
colnames(heart4)[3] <- "SAMPID" 
colnames(liver4)[3] <- "SAMPID" 

# Add sample ID col to the long-format datatframes
heart4$SUBJID <- str_sub(heart4$SAMPID, start=1L, end=10L)
liver4$SUBJID <- str_sub(liver4$SAMPID, start=1L, end=10L)

head(heart4)
head(liver4)

# There should be 117 * 296 rows = 17316 
dim(heart4) # 17316 5
dim(liver4) # 17316 5

length(unique(heart4$gene)) # 117
any(is.na(heart4))==TRUE # FALSE

# Add column with the expression batch ID (SMGEBTCH) and the type of genotype
# or expression batch (SMGEBTCHT)
heart4$SMGEBTCH <- manifest$SMGEBTCH[match(heart4$SAMPID, manifest$SAMPID)]
liver4$SMGEBTCH <- manifest$SMGEBTCH[match(liver4$SAMPID, manifest$SAMPID)]

heart4$SMGEBTCHT <- manifest$SMGEBTCHT[match(heart4$SAMPID, manifest$SAMPID)]
liver4$SMGEBTCHT <- manifest$SMGEBTCHT[match(liver4$SAMPID, manifest$SAMPID)]

heart4$SMRIN <- manifest$SMRIN[match(heart4$SAMPID, manifest$SAMPID)]
liver4$SMRIN <- manifest$SMRIN[match(liver4$SAMPID, manifest$SAMPID)]

heart4$SMTSISCH <- manifest$SMTSISCH[match(heart4$SAMPID, manifest$SAMPID)]
liver4$SMTSISCH <- manifest$SMTSISCH[match(liver4$SAMPID, manifest$SAMPID)]

head(heart4)
head(liver4)

# Try adding an index to each dataframe to combine them
heart4 <- heart4 %>% mutate(id=row_number())
liver4 <- liver4 %>% mutate(id=row_number())

# Convert subject_ID column to factor
heart4$id <- as.factor(heart4$id)
liver4$id <- as.factor(liver4$id)

# There should be 117 genes * 148 samples = 17316 rows
dim(heart4) # 17316 10
dim(liver4)

# Dataframe with only samples from people who donated both a heart and liver: organs
# e.g. paired samples
qnorm_voom_SLC25_biomarkers_manifest <- bind_rows(heart4, liver4) 
#qnorm_voom_SLC25_biomarkers_manifest <- left_join(heart4, liver4, by="SUBJID") 
length(unique(qnorm_voom_SLC25_biomarkers_manifest$SUBJID)) # 148 individuals

# Check if any genes are missing
vec <- qnorm_voom_SLC25_biomarkers_manifest$gene
vec1 <- unique(vec)
vec1 %in% genes_to_subset
vec1[118] # ID; there is an extra row from binding... 

# There should be 53 SLC25 + 117 biomarkers * 148 individuals * 2 samples per individual = 34632 rows
# Rename the columns to be more readable
cols <- c("gene", "organ", "sample_ID", "qnorm_voom_logCPM", "subject_ID", "batch_ID", "sequencing_machine", "RIN", "ischemic_time")
colnames(qnorm_voom_SLC25_biomarkers_manifest) <- cols
head(qnorm_voom_SLC25_biomarkers_manifest)
tail(qnorm_voom_SLC25_biomarkers_manifest)
dim(qnorm_voom_SLC25_biomarkers_manifest) #34632 10 

# Check for duplicated rows
# There should be 53 SLC25 * 117 biomarkers * 148 * 2
dim(unique(qnorm_voom_SLC25_biomarkers_manifest)) #34632 10 

write.table(qnorm_voom_SLC25_biomarkers_manifest, "/scratch/mjpete11/mitochondrial_paralogs/linear_models/data/data/qnorm_voom_SLC25_biomarkers_manifest.csv", sep=",")
