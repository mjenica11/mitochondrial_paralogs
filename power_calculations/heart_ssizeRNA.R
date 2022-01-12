library(ssizeRNA)
library(edgeR)
library(stringr)

# Read in GTEx metadata
meta <- read.delim("/scratch/mjpete11/power_calculations/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt",header = TRUE, sep = "\t")

# Keep only heart samples
heart <- meta[meta$SMTS == "Heart",]

# Filter rows with NA in the RIN column
test <- heart[!is.na(heart$SMRIN),]

# Drop samples with RIN < 8.0
heart <- heart[heart$SMRIN > 6.0,]

# Read in GTEx count data
counts <- read.delim("/scratch/mjpete11/power_calculations/quantification_gencode.counts.txt", header = TRUE, sep = "\t")

# Replace . to - in colnames
colnames(counts) <- str_replace_all(colnames(counts), pattern="\\.","-")

# Subset counts to only samples present in metadata
test <- counts[colnames(counts)[!names(counts) %in% heart$SAMPID]]

# Convert df to matrix
counts <- data.matrix(counts)

# Filter zero counts genes
counts <- counts[rowSums(counts) > 0,]

# Average read count for each gene
mu <- apply(counts, 1, mean)

# Dispersion for each gene
d <- DGEList(counts)
d <- calcNormFactors(d)
d <- estimateCommonDisp(d)
d <- estimateTagwiseDisp(d)
disp <- d$tagwise.dispersion

# Estimate sample size 
set.seed(2022)
pdf("heart_power_calc.pdf")
size <- ssizeRNA_vary(nGenes = 50, pi0 = 0.8, m = 1000, mu = mu, disp = disp,
					  fc = 2, fdr = 0.05, power = 0.8, maxN = 150, replace = FALSE)
dev.off()
