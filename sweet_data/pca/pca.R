# PCA plots of unadjusted, quantile normalize, and combat normalized sweet count data

library(tidyr)
library(data.table)
library(dplyr)
library(stats)
library(ggplot2)
library(ggfortify)
library(tidyverse)
library(stringr)
library(zoo)
library(readxl)
library(janitor)
library(viridis)

# Path
path <- getwd()
path

## Read in filtered not not normalized data 
# Read in clinical data 
nf_counts <- read.table("/scratch/mjpete11/mitochondrial_paralogs/sweet_data/count_matrices/filtered_non_failing_controls.csv", sep=",")
icm_counts <- read.table("/scratch/mjpete11/mitochondrial_paralogs/sweet_data/count_matrices/filtered_ischemic_cardiomyopathy.csv", sep=",")
dcm_counts <- read.table("/scratch/mjpete11/mitochondrial_paralogs/sweet_data/count_matrices/filtered_dilated_cardiomyopathy.csv", sep=",")

# Combine into one dataframe
nf_counts[1:5,1:5]
icm_counts[1:5,1:5]
dcm_counts[1:5,1:5]

# Dimensions
dim(nf_counts) # 24210    16 
dim(icm_counts) # 24688    15  
dim(dcm_counts) # 25571    38 --> 63 samples total

# Make the first row into the colnames
nf_counts <- nf_counts %>% row_to_names(row_number=1) 
nf_counts[1:5,1:5]
colnames(nf_counts)
colnames(nf_counts)[1] <- "index"

icm_counts <- icm_counts %>% row_to_names(row_number=1) 
icm_counts[1:5,1:5]
colnames(icm_counts)
colnames(icm_counts)[1] <- "index"

dcm_counts <- dcm_counts %>% row_to_names(row_number=1) 
dcm_counts[1:5,1:5]
colnames(dcm_counts)
colnames(dcm_counts)[1] <- "index"

# Combine into one dataframe
counts <- list(nf_counts, icm_counts, dcm_counts) %>% reduce(inner_join, by=c("Hugo_ID"))
head(counts)
dim(counts) #  23346 65  --> 63 samples
colnames(counts)
head(counts$Hugo_ID) 

# Subset to the top 500 most highly expressed genes
dim(counts) # 23346 67
counts[1:5, 1:5]
class(as.data.frame(counts))
apply(counts, 2, class)
counts_numeric <- apply(counts[,3:ncol(counts)], 2, as.numeric)
apply(counts_numeric, 2, class)
class(counts_numeric)
counts_numeric[1:5,1:5]
min(colMeans(counts_numeric)) # 828
max(colMeans(counts_numeric)) # 1464 
#select = order(rowMeans(counts_numeric[,3:ncol(counts_numeric)]), decreasing=TRUE)[1:500]
#highexprgenes_counts <- counts_numeric[,3:ncol(counts_numeric)][select,]
select = order(rowMeans(counts_numeric), decreasing=TRUE)[1:500]
highexprgenes_counts <- counts_numeric[select,]

#dim(highexprgenes_counts) # 500 64 
dim(highexprgenes_counts) # 500 63 
apply(highexprgenes_counts, 2, class)
highexprgenes_counts[1:5, 1:5]

# Check the dimensions
dim(highexprgenes_counts) # 500 63 

tail(colnames(highexprgenes_counts))
head(colnames(highexprgenes_counts))
highexprgenes_counts[complete.cases(highexprgenes_counts), ]
head(highexprgenes_counts)
# append column names back on
colnames(highexprgenes_counts) <- colnames(counts)[3:ncol(counts)]
counts1 <- highexprgenes_counts

# Drop the extraneous columns
counts1 <- data.frame(counts1) %>% select(contains("SAMN"))
colnames(counts1)
dim(counts1) # 500 64
counts1[1:5,1:5]

# Check the class of the columns
apply(counts1, 2, class)
#counts1 <- apply(counts1[,3:ncol(counts1)], 2, as.numeric)
counts1 <- apply(counts1, 2, as.numeric)
head(apply(counts1, 2, class))
counts1[1:5,1:5]
any(is.na(counts1))==TRUE # FALSE 
dim(counts1) # 24147    64

# Read in sample attributes files
meta <- read.csv("/scratch/mjpete11/mitochondrial_paralogs/sweet_data/count_matrices/subset_metadata.csv", sep=",")
dim(meta) # 64 30; two samples are missing from the counts1 data

# Read in the second metadata file with the RIN scores
#metadata <- read_excel("/scratch/mjpete11/mitochondrial_paralogs/sweet_data/count_matrices/sweet_metadata2.xlsx")
#dim(metadata) # 65 7
#
## Drop the 27th sample; this one couldn't be processed 
#class(metadata)
#metadata <- data.frame(metadata)
#metadata <- row_to_names(metadata, row_number=1)
#dim(metadata) # 64 7
#head(metadata)
##metadata <- metadata[-c(27),]
##dim(metadata) # 63 7
#
#metadata[1:5,1:7]
#colnames(metadata)
#metadata$RIN <- as.numeric(metadata[[7]])
#metadata$RIN
#
## Add the entire metadata for the samples in the corresponding umap projections
#meta <- manifest[manifest$'BioSample' %in% colnames(counts1),]
#meta <- meta[, c(1,2,7,14,25,28)]
#meta <- meta[-c(27),]
head(meta)
dim(meta) # 64 6

# Check class factor
#apply(meta, 2, class)
#class(meta)
#
## Factor response variables
#meta$disease <- as.factor(meta$disease)
#meta$sex <- as.factor(meta$sex)
#meta$AGE <- as.numeric(meta$AGE)
#
## Check class factor
#apply(meta, 2, class)
#
## Add the RIN
#nrow(meta) # 64
#nrow(metadata) # 65
#nrow(manifest) # 64
#length(colnames(counts1))
#
## Which sample is missing; This sample couldn't be processed by Salmon
#setdiff(manifest$BioSample, meta$BioSample) # "SAMN09484141"
#
## Add the RIN score
#meta$RIN <- metadata$RIN
#meta$RIN
#
## Order the samples in the count data to be in the same order
## as the meta metadata so they are labelled correctly
#head(meta)
idx <- match(meta$BioSample, colnames(counts1))
counts12 <- as.data.table(counts1)
ordered_counts1 <- counts12[,..idx]
all(meta$BioSample==colnames(ordered_counts1)) # TRUE

# PCA on unnormalized data
# Set seed
set.seed(1993)
ordered_counts1[1:5,1:5]
#pca_obj <- prcomp(t(ordered_counts1[,3:ncol(ordered_counts1)]))
#pca_obj <- prcomp(t(ordered_counts1))
pca_obj <- prcomp(ordered_counts1) # skipping the transposition
class(pca_obj)
names(pca_obj)
head(pca_obj$x)
head(ordered_counts1)
head(pca_obj$x[,1])
head(pca_obj$x[,2])

# Save 1st and 2nd dimensions as objects
#pc1 <- pca_obj$x[,1]
#pc2 <- pca_obj$x[,2]
#colnames(pc1)==colnames(pc2)

# Regress out PC1 and PC2
#pc_adjusted_lm_obj <- lm(pca_obj$x~pc1+pc2)
#dim(pc_adjusted_lm_obj)
#head(pc_adjusted_lm_obj)
#class(pc_adjusted_lm_obj)
#names(pc_adjusted_lm_obj)
#hist(pc_adjusted_lm_obj$effects, breaks=100)
#pc_adjusted <- as.data.frame(pc_adjusted_lm_obj$effects)
#class(pc_adjusted)
#dim(pc_adjusted) # 63 63

# Check the data was successfully adjusted
#pc_adjusted[1:5,1:5]
#pca_obj$x[1:5,1:5]

# Call the object pc_adjusted even if I didn't adjust
# so I don't have to rename all the following objects
# The object is transposed (gene x sample) in order
# to facilitate merging with the metadata
# Specifically, the columns need to be the rows 
pc_adjusted <- as.data.frame(t(pca_obj$x))
class(pc_adjusted)
dim(pc_adjusted) # 63 500 

# Add the rownames as a column
pc_adjusted$BioSample <- meta$BioSample
head(pc_adjusted)
class(pc_adjusted)
tail(pc_adjusted)
pc_adjusted[1:5,1:4]
dim(pc_adjusted) # 63 501 
str(pc_adjusted)

# Merge the pc_adjusted  projections with the metadata
dim(meta) # 63 7
head(pc_adjusted)

# Make a copy of the meta object without the samples missing in the count data
meta_sub <- meta[meta$BioSample %in% pc_adjusted$BioSample,]
dim(meta_sub) # 63 7
all(meta_sub$BioSample==pc_adjusted$BioSample) # TRUE
idx2 <- match(meta_sub$BioSample, pc_adjusted$BioSample)
pc_adjusted2 <- as.data.table(pc_adjusted)
ordered_pc_adjusted <- pc_adjusted2[,..idx2]
ordered_pc_adjusted[1:5,1:5]
colnames(ordered_pc_adjusted)
dim(ordered_pc_adjusted) # 63 63
# Add the sample names on as a column so you can merge
ordered_pc_adjusted$BioSample <- meta_sub$BioSample
ordered_pc_adjusted$BioSample
all(meta_sub$BioSample==ordered_pc_adjusted$BioSample) # TRUE
# merge
pc_adjusted <- merge(ordered_pc_adjusted, meta_sub, by="BioSample")
head(pc_adjusted)
dim(pc_adjusted) # 63 70
colnames(pc_adjusted)
class(pc_adjusted)
pc_adjusted <- data.frame(pc_adjusted)
class(pc_adjusted)
head(pc_adjusted)

# Two samples had projections that were dramatically different
WeirdSamples<-c("SAMN09484136","SAMN09484137")

pc_adjusted%>%
  filter(!BioSample%in%WeirdSamples)%>%
  #filter(PC1>0)%>%
  ggplot(aes(PC1,PC2,color = disease))+
  geom_point()

## Save the pca data as a df
#class(pc_adjusted)
#pca_df <- as.data.frame(pc_adjusted)
##rownames(pca_obj$x)
#rownames(pc_adjusted)
#rownames(pca_df)
#colnames(pc_adjusted)
#colnames(pca_df)
#
## Add the rownames as a column 
#pca_df$BioSample <- rownames(pca_obj$x)
#class(pca_df)
#tail(pca_df)
#head(pca_df)
#
## Merge the pca_dfections with the metadata
#dim(pca_df) # 63 64
#dim(meta) # 63 7
## Make a copy of the meta object without the samples missing in the count data
#meta_sub <- meta[meta$BioSample %in% pca_df$BioSample,]
#dim(meta_sub) # 62 7
#all(meta_sub$BioSample==pca_df$BioSample) # TRUE 
#idx2 <- match(meta_sub$BioSample, pca_df$BioSample)
#pca_df2 <- as.data.table(pca_df)
#ordered_pca_df <- pca_df2[,..idx2]
#ordered_pca_df[1:5,1:5]
#colnames(ordered_pca_df)
#dim(ordered_pca_df) # 63 63
## Add the sample names on as a column so you can merge
#ordered_pca_df$BioSample <- meta_sub$BioSample
#ordered_pca_df$BioSample
#all(meta_sub$BioSample==ordered_pca_df$BioSample) # TRUE
## merge
#pca_df <- merge(ordered_pca_df, meta_sub, by="BioSample")
#head(pca_df)
#
# Plot PCA colored by RIN 
#p1 <- ggplot(pca_df, aes(PC1,PC2,color=RIN)) + geom_point(size=1, show.legend=TRUE) + scale_color_viridis()

# Print plot
#pdf(file=paste0(path,"/batch_normalized_PCA.pdf"))
#pdf(file=paste0(path,"/RIN_filtered_PCA.pdf"))
#p1
#dev.off()

# Plot PCA colored by sex 
p2 <- ggplot(pc_adjusted, aes(V1,V2,color=sex)) + geom_point(size=1, show.legend=TRUE) + geom_jitter() + labs(x="PC1", y="PC2")
# Print plot
#png(file=paste0(path,"/sex_pc_adjusted_PCA.png"))
png(file=paste0(path,"/sex_filtered_PCA.png"))
p2
dev.off()

# Plot PCA colored by disease 
p3 <- ggplot(pc_adjusted, aes(V1,V2,color=disease)) + geom_point(size=1, show.legend=TRUE) + geom_jitter() + labs(x="PC1", y="PC2")
# Print plot
#png(file=paste0(path,"/disease_pc_adjusted_PCA.png"))
png(file=paste0(path,"/disease_filtered_PCA.png"))
p3
dev.off()

# Plot PCA colored by age 
p4 <- ggplot(pc_adjusted, aes(V1,V2,color=AGE)) + geom_point(size=1, show.legend=TRUE) + scale_color_viridis() + geom_jitter() + labs(x="PC1", y="PC2")
# Print plot
#png(file=paste0(path,"/age_pc_adjusted_PCA.png"))
png(file=paste0(path,"/age_filtered_PCA.png"))
p4
dev.off()


# Calculate percent variance
print(head(((pca_obj$sdev^2) / (sum(pca_obj$sdev^2)))*100))

