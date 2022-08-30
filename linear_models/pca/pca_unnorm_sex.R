# PCA plots of unadjusted, quantile normalize, and combat normalized GTEx count data

library(tidyr)
library(data.table)
library(dplyr)
library(stats)
library(ggplot2)
library(ggfortify)
library(tidyverse)
library(stringr)

# Read in quantile normalized GTEx data
counts <- fread("/scratch/mjpete11/linear_models/data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct", sep="\t")
#counts <- counts[,1:5000]

# Read in metadata with batch variables
manifest <- read.csv("/scratch/mjpete11/linear_models/data/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt", header=TRUE, sep = "\t")

# Add sex and age to manifest
# Read in separate file
meta <- read.csv("/scratch/mjpete11/linear_models/data/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt", header=TRUE, sep = "\t")

# Add sample ID col to meta
meta$SAMPID <- manifest$SAMPID[sapply(meta$SUBJID, function(x) match(x, substr(manifest$SAMPID, 1, nchar(x))))]

# Get list of female IDs
fems <- meta$SUBJID[which(meta$SEX==2)]

# Make list of sex of each individual
sex <- with(meta['SUBJID'], ifelse(SUBJID %in% fems, "Female", "Male"))

# Add column containing sex
meta <- cbind(meta, "Sex"=sex)

# Add column containing age (only decade intervals are publically accessible)
#meta$Age <- meta$AGE[match(meta$Individual_ID, meta$SUBJID)]

# Bind sex and age columns to manifest by matching sample IDs
manifest <- right_join(meta, manifest, by=c("SAMPID"))

# Remove all rows with NA in either columns
#manifest <- manifest %>% drop_na(SAMPID, SMTSD, SMTSISCH, SMGEBTCH, SMRIN, SEX, AGE) 
manifest <- manifest %>% drop_na(SAMPID, SMTSD) 

# Does the number of samples in the manifest == number of samples in counts?
nrow(manifest)==ncol(counts[3:ncol(counts)])

# Median filter: Function to filter rows that do not have a median >= condition
median_filter <- function(DF, thresh){
		DF <- DF[which(apply(DF[,-c(1:2)],1,median) > thresh), ]
		return(DF)
}

# Filter rows with a median of < 10 counts
counts <- median_filter(DF=counts, thresh=10)

# Drop samples from manifest that are missing in the count data
#rownames(manifest) <- manifest$SAMPID
manifest <- manifest[(manifest$SAMPID %in% colnames(counts)[3:ncol(counts)]),]
samples <- manifest$SAMPID

# Only keep samples that are present in the manifest 
counts <- counts[,..samples]

# TEMP
#counts <- counts[,3:500]

# Subset heart left ventricle and liver samples into two separate dfs
heart <- manifest[manifest$SMTSD %in% "Heart - Left Ventricle", ]
liver <- manifest[manifest$SMTSD %in% "Liver", ]

# Subset the SLC25 gene count df by the sample IDs that match the IDs in file3 df
heart2 <- counts %>% select(contains(heart$SAMPID))
liver2 <- counts %>% select(contains(liver$SAMPID))

# Append the gene ensemble ID and common name
heart3 <- cbind('gene'=counts$'Description', heart2)
liver3 <- cbind('gene'=counts$'Description', liver2)

heart3 <- cbind('gene'=counts[,1], heart2)
liver3 <- cbind('gene'=counts[,1], liver2)
colnames(heart3)[1] <- 'gene' 
colnames(liver3)[1] <- 'gene' 

# Add a columns with the organ in each df
heart3$organ <- 'heart' 
liver3$organ <- 'liver' 

# Reshape dataframe so it can be converted to a design matrix object
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

heart4$SEX <- manifest$SEX[heart4$SAMPID %in% manifest$SAMPID]
liver4$SEX <- manifest$SEX[liver4$SAMPID %in% manifest$SAMPID]

heart4$AGE <- manifest$AGE[heart4$SAMPID %in% manifest$SAMPID]
liver4$AGE <- manifest$AGE[liver4$SAMPID %in% manifest$SAMPID]

# Combine into one df
organs <- rbind(heart4, liver4)

# Factor response variable
organs$organ <- as.factor(organs$organ)

# PCA on unnormalized data
#pca_obj <- prcomp(t(counts))
pca_obj <- prcomp(t(counts))

# PCA 
df <- data.frame(pca_obj$x[,1:2])
df$SAMPID <- row.names(pca_obj$x)
df_2 <- merge(organs %>% select(organ,SAMPID,SMGEBTCH, SMTSISCH, SMRIN, SEX, AGE) %>% unique(), df)

#p1 <- ggplot(df_2, aes(PC1,PC2,color=SMTSISCH)) + geom_point(show.legend=FALSE)
p1 <- ggplot(df_2, aes(PC1,PC2,color=SEX)) + geom_point(show.legend=TRUE)

# Calculate percent variance
print(head(((pca_obj$sdev^2) / (sum(pca_obj$sdev^2)))*100))

# Print plot
pdf(file="/scratch/mjpete11/linear_models/results2/sex_unnormalized_PCA.pdf")
p1
dev.off()
