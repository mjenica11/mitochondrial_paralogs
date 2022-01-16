library(ssizeRNA)
library(edgeR)
library(stringr)
library(data.table)
library(stringr)

# Read in GTEx metadata
meta <- read.delim("/scratch/mjpete11/power_calculations/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt",header = TRUE, sep = "\t")

# Read in file with sex phenotypes
sex_meta <- read.delim("/scratch/mjpete11/power_calculations/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt",header = TRUE, sep = "\t")

# Get list of female IDs
fems <- sex_meta$SUBJID[which(sex_meta$SEX==2)]

# Add columns with the invidual IDs and the sex of each sample 
meta[["Individual_ID"]] <- str_extract(meta$SAMPID, "GTEX-[0-9A-Z]+")
sex <- with(meta['Individual_ID'], ifelse(Individual_ID %in% fems, "Female", "Male"))
meta[['sex']] <- sex

# Dispersion of heart counts only
meta2 <- meta[which(meta$SMTS == "Heart"),]

# Filter rows with NA in the RIN column
meta2 <- meta2[!is.na(meta2$SMRIN),]

# Drop samples with RIN < 8.0
meta2 <- meta2[which(meta2$SMRIN > 8.0),]

# Read in GTEx count data
#counts <- data.frame(fread("/scratch/mjpete11/power_calculations/GTEx_v8_counts.gct"))

# Replace . to - in colnames
#colnames(counts) <- str_replace_all(colnames(counts), pattern="\\.","-")

# Subset count matrix to only the heart samples
#subcounts  <- counts[colnames(counts) %in% meta2$SAMPID]

# Add gene name column
#subcounts <- cbind(gene=counts$Name, subcounts)

# Write subsetted counts to disk
#write.csv(subcounts, "/scratch/mjpete11/power_calculations/heart_counts.csv")

# Read in the heart sample count data only
subcounts <- data.frame(fread("/scratch/mjpete11/power_calculations/heart_counts.csv"))

# Replace . to - in colnames
colnames(subcounts) <- str_replace_all(colnames(subcounts), pattern="\\.","-")

# Convert df to matrix
subcounts <- data.matrix(subcounts)

# Filter zero counts genes
subcounts <- subcounts[rowSums(subcounts) > 0,]

# Average read count for each gene
mu <- apply(subcounts, 1, mean)

# Drop samples from metadata that are missing count data
select_samples <- colnames(subcounts)[colnames(subcounts) %in% meta2$SAMPID]
meta3 <- meta2[meta2$SAMPID %in% select_samples,]

# Drop first two cols from count matrix
subcounts <- subcounts[,3:ncol(subcounts)]

# Dispersion for each gene
d <- DGEList(subcounts)
d <- calcNormFactors(d)
d <- estimateCommonDisp(d)
d <- estimateTagwiseDisp(d)
disp <- d$tagwise.dispersion

# Estimate sample size 
set.seed(2022)
pdf("heart_power_calc.pdf")
size <- ssizeRNA_vary(nGenes = 50, pi0 = 0.2, m = 50, mu = mu, disp = disp,
					  fc = 2, fdr = 0.05, power = 0.8, maxN = 15, replace = FALSE)
dev.off()

# Differential expression between male and female samples
sex_groups <- factor(meta3$sex)
d2 <- DGEList(subcounts, group=sex_groups)
mat <- model.matrix(~group, data=d2[['samples']])
v <- voom(d2, mat)
fit <- lmFit(v, mat)
fit <- eBayes(fit)

# Apply FDR correction and then drop p-vals with >0.05
res <- topTable(fit, adjust.method="fdr", number=nrow(fit))
res2 <- subset(res, adj.P.Val < 0.05)

# Number of significantly DE genes after multiple test correction
nrow(res2) # 342

# Total number of genes tested with expression > 0 in at least 1 sample
nrow(subcounts) # 51,464

# Proportion of genes expected to be DE between healthy fe/male heart tissue
342/51464 #0.0067

# Dispersion for each gene in sex-stratified heart samples
d3 <- calcNormFactors(d2)
d3 <- estimateCommonDisp(d3)
d3 <- estimateTagwiseDisp(d3)
disp2 <- d3$tagwise.dispersion

# Average read count for each gene
mu2 <- apply(subcounts, 1, mean)

# Estimate sample size 
set.seed(2022)
pdf("sex_heart_power_calc.pdf")
size2 <- ssizeRNA_vary(nGenes=51500, pi0=0.007, m=50, mu=mu2, disp=disp2,
	   				   fc=2, fdr=0.05, power=0.8, maxN=15, replace=TRUE)
dev.off()
