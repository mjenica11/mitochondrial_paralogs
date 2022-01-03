# Differential expression of simulated SLC25 transcripts

# Load packages
library(limma)
library(edgeR)
library(data.table)

# Constants
GENE_COUNTS <- "~/Downloads/gene_count_matrix.tsv"

# Read in counts
gene_counts <- read.csv(GENE_COUNTS, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
head(gene_counts)

#_______________________________________________________________________________________ 
# Make dataframe
#_______________________________________________________________________________________ 
# Subset the remaining SLC25 genes from the count matrix
SLC25A4 <- as.numeric(gene_counts[9160,])
SLC25A5 <- as.numeric(gene_counts[76,]) 
SLC25A6 <- as.numeric(gene_counts[12103,]) 
SLC25A31 <- as.numeric(gene_counts[9121,]) 
SLC25A43 <- as.numeric(gene_counts[1376,]) 
SLC25A24 <- as.numeric(gene_counts[1635,]) 
SLC25A25 <- as.numeric(gene_counts[8803,]) 
SLC25A23 <- as.numeric(gene_counts[5551,]) 
SLC25A41 <- as.numeric(gene_counts[14313,]) 
SLC25A16 <- as.numeric(gene_counts[5219,]) 
SLC25A42 <- as.numeric(gene_counts[14289,]) 
SLC25A44 <- as.numeric(gene_counts[10202,]) 
SLC25A49 <- as.numeric(gene_counts[7259,]) 
SLC25A50 <- as.numeric(gene_counts[3651,]) 
SLC25A3 <- as.numeric(gene_counts[1289,]) 
SLC25A7 <- as.numeric(gene_counts[3599,]) 
SLC25A8 <- as.numeric(gene_counts[13396,]) 
SLC25A9 <- as.numeric(gene_counts[13395,]) 
SLC25A27 <- as.numeric(gene_counts[9350,]) 
SLC25A14 <- as.numeric(gene_counts[2580,]) 
SLC25A30 <- as.numeric(gene_counts[13136,]) 
SLC25A34 <- as.numeric(gene_counts[10406,]) 
SLC25A35 <- as.numeric(gene_counts[5519,]) 
SLC25A10 <- as.numeric(gene_counts[14635,]) 
SLC25A11 <- as.numeric(gene_counts[3484,]) 
SLC25A53 <- as.numeric(gene_counts[31351,]) 
SLC25A51 <- as.numeric(gene_counts[5193,]) 
SLC25A29 <- as.numeric(gene_counts[16292,]) 
SLC25A47 <- as.numeric(gene_counts[7698,]) 
SLC25A45 <- as.numeric(gene_counts[10363,]) 
SLC25A48 <- as.numeric(gene_counts[8507,]) 
SLC25A20 <- as.numeric(gene_counts[13879,]) 
SLC25A15 <- as.numeric(gene_counts[2649,]) 
SLC25A32 <- as.numeric(gene_counts[11117,]) 
SLC25A17 <- as.numeric(gene_counts[2252,]) 
SLC25A46 <- as.numeric(gene_counts[10931,]) 
SLC25A33 <- as.numeric(gene_counts[12650,]) 
SLC25A36 <- as.numeric(gene_counts[4181,]) 
SLC25A1 <- as.numeric(gene_counts[2147,]) 
SLC25A12 <- as.numeric(gene_counts[4414,]) 
SLC25A13 <- as.numeric(gene_counts[66,]) 
SLC25A18 <- as.numeric(gene_counts[14602,]) 
SLC25A22 <- as.numeric(gene_counts[13734,]) 
SLC25A21 <- as.numeric(gene_counts[14630,]) 
SLC25A26 <- as.numeric(gene_counts[8368,]) 
SLC25A28 <- as.numeric(gene_counts[9557,]) 
SLC25A37 <- as.numeric(gene_counts[8714,]) 
SLC25A38 <- as.numeric(gene_counts[8356,]) 
SLC25A39 <- as.numeric(gene_counts[326,]) 
SLC25A40 <- as.numeric(gene_counts[1281,]) 
SLC25A19 <- as.numeric(gene_counts[5524,]) 

# Add the remaining SLC25 genes to the metadata 
# SLC24A2 and SLC25A52 are missing from the count table 
lst <- list(SLC25A4,SLC25A5,SLC25A6,SLC25A31,SLC25A43,SLC25A24,SLC25A25,
			SLC25A23,SLC25A41,SLC25A16,SLC25A42,SLC25A44,SLC25A49,SLC25A50,
			SLC25A3,SLC25A7,SLC25A8,SLC25A9,SLC25A27,SLC25A14,SLC25A30,
			SLC25A34,SLC25A35,SLC25A10,SLC25A11,SLC25A53,SLC25A51,SLC25A29,
			SLC25A47,SLC25A45,SLC25A48,SLC25A20,SLC25A15,SLC25A32,SLC25A17,
			SLC25A46,SLC25A33,SLC25A36,SLC25A1,SLC25A12,SLC25A13,SLC25A18,
			SLC25A22,SLC25A21,SLC25A26,SLC25A28,SLC25A37,SLC25A38,SLC25A39,
			SLC25A40,SLC25A19)
dat <- t(data.frame(matrix(unlist(lst), nrow=length(lst), byrow=TRUE)))
names <- c("SCL25A4","SLC25A5","SLC25A6","SLC25A31","SLC25A43","SLC25A24",
		   "SLC25A25","SLC25A23","SLC25A41","SLC25A16","SLC25A42","SLC25A44",
		   "SLC25A49","SLC25A50","SLC25A3","SLC25A7","SLC25A8","SLC25A9",
		   "SLC25A27","SLC25A14","SLC25A30","SLC25A34","SLC25A35","SLC25A10",
		   "SLC25A11","SLC25A53","SLC25A51","SLC25A29","SLC25A47","SLC25A45",
		   "SLC25A48","SLC25A20","SLC25A15","SLC25A32","SLC25A17","SLC25A46",
		   "SLC25A33","SLC25A36", "SLC25A1", "SLC25A12","SLC25A13","SLC25A18",
		   "SLC25A22","SLC25A21","SLC25A26","SLC25A28","SLC25A37","SLC25A38",
		   "SLC25A39","SLC25A40","SLC25A19")
colnames(dat) <- names

# Add sample IDs and condition to metadata
samples <- sprintf("sample_%02d", 01:20)
condition <- c(replicate(10, "HEALTHY"), replicate(10, "SICK"))

# Combine dataframes and drop individual and group columns
dat <- cbind(dat, data.frame(samples = samples))
dat <- cbind(dat, data.frame(condition = condition))

# Make rownames the sample names
row.names(dat) <- samples

# Write table of counts, sample ID, and condition
#write.csv(dat, "~/Downloads/SLC25_counts.csv")

#_______________________________________________________________________________________ 
# Differential expression of all SLC25 genes between conditions 
#_______________________________________________________________________________________ 
# Design matrix for all SLC25 members
mat <- model.matrix(~0 + condition, data=dat)

# Apply function to each gene to determine the p-values, betas, etc.
fit_function <- function(dat_object) {
				res <- eBayes(lmFit(dat_object, mat))
				return(res)
}
res <- apply(dat[,1:51],2,fit_function)

# Make table of p-values, standard deviation, and betas for each gene 
make_tables <- function(x){
		tab <- data.frame(group1_pval=x[["p.value"]][,1],
    					  group2_pval=x[["p.value"]][,2],
						  stdev=x[["sigma"]],
						  group1_beta=x[["coefficients"]][,1],
						  group2_beta=x[["coefficients"]][,2])
		return(tab)
}
table_lst <- lapply(res, make_tables)

# Convert list of dataframes into one dataframe
dat2 <- rbindlist(table_lst) 

# Add column with gene names
row.names(dat2) <- names

# Write table of p-values, standard deviation, and betas for each gene
#write.csv(dat3, "~/Downloads/p_values_and_betas.csv")

# Test for differential expression
res2 <- t(sapply(res, decideTests,method="separate", adjust.method="BH", p.value=0.05, lfc=2))

# Make a venn diagram of the significantly DE genes between conditions
obj <- vennCounts(res2)
print(obj)
mfrow.old <- par()$mfrow
par(mfrow=c(1,2))
vennDiagram(obj)
vennDiagram(obj, 
            include=c("up", "down"),
            counts.col=c("red", "blue"),
            circle.col=c("red", "blue", "green3"))
par(mfrow=mfrow.old)
