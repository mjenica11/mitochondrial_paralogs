# Extract the number of reads that mapped to the SLC25 genes
dat <- read.csv("featureCounts.count", sep="\t", header=TRUE)

# Genes to extract
SLC25_genes  <- c("gene-SLC25A33", "gene-SLC25A34", "gene-SLC25A24", "gene-SLC25A44", "gene-SLC25A12",
				  "gene-SLC25A38", "gene-SLC25A20", "gene-SLC25A26", "gene-SLC25A36", "gene-SLC25A31",
				  "gene-SLC25A7", "gene-SLC25A4", "gene-SLC25A46", "gene-SLC25A48", "gene-SLC25A2",
				  "gene-SLC25A49", "gene-SLC25A27", "gene-SLC25A40", "gene-SLC25A13", "gene-SLC25A37",
				  "gene-SLC25A32", "gene-SLC25A51", "gene-SLC25A25", "gene-SLC25A16", "gene-SLC25A28",
				  "gene-SLC25A50", "gene-SLC25A45", "gene-SLC25A8", "gene-SLC25A9", "gene-SLC25A22",
				  "gene-SLC25A3", "gene-SLC25A15", "gene-SLC25A30", "gene-SLC25A21", "gene-SLC25A29",
				  "gene-SLC25A47", "gene-SLC25A11", "gene-SLC25A35", "gene-SLC25A39", "gene-SLC25A19",
				  "gene-SLC25A10", "gene-SLC25A52", "gene-SLC25A41", "gene-SLC25A23", "gene-SLC25A42",
				  "gene-SLC25A18", "gene-SLC25A1", "gene-SLC25A17", "gene-SLC25A6", "gene-SLC25A53",
				  "gene-SLC25A43", "gene-SLC25A5", "gene-SLC25A14")

# Subset SLC25 genes
dat2 <- subset(dat, dat$Geneid %in% SLC25_genes)
nrow(dat2) # 48

# Drop unnecessary columns
dat2 <- dat2[,c(1,6,7)]

# Change the rownames
rownames(dat2) <- seq(1:48)

# What is the index of the gene the reads were simulated from?
which(dat2$Geneid %in% "gene-SLC25A4") # 11

# Add a column with the number of simulated reads
simulated <- rep(0, 48)
simulated[11] <- 2000
dat2$simulated_reads <- simulated 

# Write to file 
write.csv(dat2, "SLC25_featureCounts.csv")
