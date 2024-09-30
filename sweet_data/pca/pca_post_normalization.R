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

# Read in quantile + voom normalized clinical data
counts <- fread("/scratch/mjpete11/mitochondrial_paralogs/sweet_data/qnorm_voom_normalize/voom_qnorm_counts.csv", sep=",")
dim(counts) # 24209 65 # --> 63 samples
counts[1:5,1:5]
colnames(counts)

# Subset to the top 500 most highly expressed genes
class(counts)
counts <- as.data.frame(counts)
apply(counts, 2, class) # everything is encoded as a character vector
counts_numeric <- apply(counts[,3:ncol(counts)], 2, as.numeric)
apply(counts_numeric, 2, class) # numeric
class(counts_numeric) # matrix array
counts_numeric[1:5,1:5]
min(colMeans(counts_numeric)) # 0.246 
max(colMeans(counts_numeric)) # 0.308 
select = order(colMeans(counts_numeric), decreasing=TRUE)[1:500]
highexprgenes_counts <- counts_numeric[select,]
apply(highexprgenes_counts, 2, class) # still numeric
highexprgenes_counts[1:5, 1:5]

# Check the dimensions
dim(highexprgenes_counts) # 500 63 

# Check if any NA
any(is.na(highexprgenes_counts)==TRUE)
highexprgenes_counts <- highexprgenes_counts[complete.cases(highexprgenes_counts), ]
head(highexprgenes_counts)
tail(highexprgenes_counts)

# Check the dimensions
dim(highexprgenes_counts) # 63 63 

# Drop the extraneous columns
counts1 <- data.frame(highexprgenes_counts) %>% select(contains("SAMN"))
colnames(counts1)
dim(counts1) # 63 64
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
ordered_counts1[1:5,1:5]
#pca_obj <- prcomp(t(ordered_counts1))
pca_obj <- prcomp(ordered_counts1)
class(pca_obj)
names(pca_obj)
head(pca_obj$x)
head(pca_obj$x[,1])
head(pca_obj$x[,2])

# Save 1st and 2nd dimensions as objects
#pc1 <- pca_obj$x[,1]
#pc2 <- pca_obj$x[,2]
#colnames(pc1)==colnames(pc2)
#
## Regress out PC1 and PC2
#pc_adjusted_lm_obj <- lm(pca_obj$x~pc1+pc2)
#dim(pc_adjusted_lm_obj)
#head(pc_adjusted_lm_obj)
#class(pc_adjusted_lm_obj)
#names(pc_adjusted_lm_obj)
#hist(pc_adjusted_lm_obj$effects, breaks=100)
#pc_adjusted <- as.data.frame(pc_adjusted_lm_obj$effects)
#class(pc_adjusted)
#dim(pc_adjusted) # 63 63
#
## Check the data was successfully adjusted
#pc_adjusted[1:5,1:5]
#pca_obj$x[1:5,1:5]
#

# Copy the unadjusted untranposed pca object
# I don't want to rename all the downstream variables
pc_adjusted <- as.data.frame(pca_obj$x)
dim(pc_adjusted) # 63 63

# Add the rownames as a column
#pc_adjusted$BioSample <- rownames(pca_obj$x)
pc_adjusted$BioSample <- meta$BioSample 
class(pc_adjusted)
tail(pc_adjusted)
pc_adjusted[1:5,1:4]
dim(pc_adjusted) # 63 64
str(pc_adjusted)
#dim(pc_adjusted$`pca_obj$x`) # 63 63
#pc_adjusted$`pca_obj$x`

# Merge the pc_adjustedections with the metadata
dim(pc_adjusted) # 63 63
dim(meta) # 63 7
head(pc_adjusted)
# Make a copy of the meta object without the samples missing in the count data
meta_sub <- meta[meta$BioSample %in% pc_adjusted$BioSample,]
dim(meta_sub) # 62 7
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

# Two samples have extremely different projections
WeirdSamples<-c("SAMN09484136","SAMN09484137")

pc_adjusted%>%
  filter(!BioSample%in%WeirdSamples)%>%
  #filter(PC1>0)%>%
  ggplot(aes(PC1,PC2,color = disease))+
  geom_point()

# Merge the pca_dfections with the metadata
dim(pc_adjusted) # 63 64
dim(meta) # 63 7
# Make a copy of the meta object without the samples missing in the count data
meta_sub <- meta[meta$BioSample %in% pc_adjusted$BioSample,]
dim(meta_sub) # 62 7
all(meta_sub$BioSample==pc_adjusted$BioSample) # TRUE 
idx2 <- match(meta_sub$BioSample, pc_adjusted$BioSample)
pca_df2 <- as.data.table(pc_adjusted)
ordered_pca_df <- pca_df2[,..idx2]
ordered_pca_df[1:5,1:5]
colnames(ordered_pca_df)
dim(ordered_pca_df) # 63 63
# Add the sample names on as a column so you can merge
ordered_pca_df$BioSample <- meta_sub$BioSample
ordered_pca_df$BioSample
all(meta_sub$BioSample==ordered_pca_df$BioSample) # TRUE
# merge
pca_df <- merge(ordered_pca_df, meta_sub, by="BioSample")
head(pca_df)

# Plot PCA colored by RIN 
#p1 <- ggplot(pca_df, aes(PC1,PC2,color=RIN)) + geom_point(size=1, show.legend=TRUE) + scale_color_viridis()

# Print plot
#pdf(file=paste0(path,"/batch_normalized_PCA.pdf"))
#pdf(file=paste0(path,"/RIN_filtered_PCA.pdf"))
#p1
#dev.off()

# Plot PCA colored by sex 
p2 <- ggplot(pc_adjusted, aes(PC1,PC2,color=sex)) + geom_point(size=1, show.legend=TRUE) + geom_jitter()
# Print plot
#png(file=paste0(path,"/sex_pc_adjusted_PCA.png"))
png(file=paste0(path,"/sex_normalized_PCA.png"))
p2
dev.off()

# Plot PCA colored by disease 
p3 <- ggplot(pc_adjusted, aes(PC1,PC2,color=disease)) + geom_point(size=1, show.legend=TRUE) + geom_jitter()
# Print plot
#png(file=paste0(path,"/disease_pc_adjusted_PCA.png"))
png(file=paste0(path,"/disease_normalized_PCA.png"))
p3
dev.off()

# Plot PCA colored by age 
p4 <- ggplot(pc_adjusted, aes(PC1,PC2,color=AGE)) + geom_point(size=1, show.legend=TRUE) + scale_color_viridis() + geom_jitter()
# Print plot
#png(file=paste0(path,"/age_pc_adjusted_PCA.png"))
png(file=paste0(path,"/age_normalized_PCA.png"))
p4
dev.off()


# Calculate percent variance
print(head(((pca_obj$sdev^2) / (sum(pca_obj$sdev^2)))*100))

