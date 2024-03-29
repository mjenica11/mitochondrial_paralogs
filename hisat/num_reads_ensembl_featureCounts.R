# Extract the number of reads that mapped to the SLC25 genes
dat <- read.csv("ensembl_featureCounts2.count", sep="\t", header=TRUE, skip=1)

# Genes to extract
SLC25_genes  <- c("SLC25A1", "SLC25A2", "SLC25A3", "SLC25A4", "SLC25A5",
				  "SLC25A6", "SLC25A7", "SLC25A8", "SLC25A9", "SLC25A10",
				  "SLC25A11", "SLC25A12", "SLC25A13", "SLC25A14", "SLC25A15",
				  "SLC25A16", "SLC25A17", "SLC25A18", "SLC25A19", "SLC25A20",
				  "SLC25A21", "SLC25A22", "SLC25A23", "SLC25A24", "SLC25A25",
				  "SLC25A26", "SLC25A27", "SLC25A28", "SLC25A29", "SLC25A30",
				  "SLC25A31", "SLC25A32", "SLC25A33", "SLC25A34", "SLC25A35",
				  "SLC25A36", "SLC25A37", "SLC25A38", "SLC25A39", "SLC25A40",
				  "SLC25A41", "SLC25A42", "SLC25A43", "SLC25A44", "SLC25A45",
				  "SLC25A46", "SLC25A47", "SLC25A48", "SLC25A49", "SLC25A50",
				  "SLC25A51", "SLC25A52", "SLC25A53")

#SLC25 <- c("ENST00000215882.10","ENST00000239451.7","ENST00000552981.6",
#		   "ENST00000281456.11","ENST00000317881.9","ENST00000381401.11",
#		   "ENST00000262999.4","ENST00000663595.2", "ENST00000314032.9",
#		   "ENST00000350690.10", "ENST00000225665.12", "ENST00000422440.7",
#		   "ENST00000265631.10", "ENST00000545805.6", "ENST00000338625.9",
#		   "ENST00000609923.6", "ENST00000435456.7", "ENST00000327451.11",
#		   "ENST00000416858.7", "ENST00000319017.5", "ENST00000331299.6",
#		   "ENST00000628067.3", "ENST00000301454.9", "ENST00000565488.6",
#		   "ENST00000373069.10", "ENST00000354883.11", "ENST00000371347.10",
#		   "ENST00000370495.6", "ENST00000359232.8", "ENST00000519676.6",
#		   "ENST00000281154.6", "ENST00000297578.9", "ENST00000302692.7",
#		   "ENST00000294454.6", "ENST00000577745.2", "ENST00000324194.12",
#		   "ENST00000519973.6", "ENST00000650617.1", "ENST00000377095.10",
#		   "ENST00000341119.10", "ENST00000321510.7", "ENST00000318596.8",
#		   "ENST00000217909.8", "ENST00000359511.5", "ENST00000398802.6",
#		   "ENST00000355943.8", "ENST00000361529.5", "ENST00000681962.1",
#		   "ENST00000373627.10", "ENST00000302503.8", "ENST00000242275.7",
#		   "ENST00000269205.7", "ENST00000594199.3")

SLC25 <- c("ENSG00000100075.10", "ENSG00000120329.7", "ENSG00000075415.15",
		   "ENSG00000151729.11", "ENSG00000005022.6", "ENSG00000169100.14",
		   "ENSG00000109424.4", "ENSG00000175567.11", "ENSG00000175564.13",
		   "ENSG00000183048.12", "ENSG00000108528.14", "ENSG00000115840.14",
		   "ENSG00000004864.14", "ENSG00000102078.16", "ENSG00000102743.15",
		   "ENSG00000122912.15", "ENSG00000100372.15", "ENSG00000182902.14",
		   "ENSG00000125454.12", "ENSG00000178537.10", "ENSG00000183032.12",
		   "ENSG00000177542.11", "ENSG00000125648.15", "ENSG00000085491.17",
		   "ENSG00000148339.13", "ENSG00000144741.19", "ENSG00000153291.16",
		   "ENSG00000155287.11", "ENSG00000197119.13", "ENSG00000174032.17",
		   "ENSG00000151475.6", "ENSG00000164933.12", "ENSG00000171612.7",
		   "ENSG00000162461.8", "ENSG00000125434.11", "ENSG00000114120.14",
		   "ENSG00000147454.14", "ENSG00000144659.13", "ENSG00000013306.16",
		   "ENSG00000075303.13", "ENSG00000181240.14", "ENSG00000181035.14",
		   "ENSG00000077713.19", "ENSG00000160785.14", "ENSG00000162241.13",
		   "ENSG00000164209.17", "ENSG00000140107.12", "ENSG00000145832.15",
		   "ENSG00000137409.20", "ENSG00000109919.10", "ENSG00000122696.14",
		   "ENSG00000141437.10", "ENSG00000269743.3")

# Subset SLC25 genes
dat2 <- subset(dat, dat$Geneid %in% SLC25)
nrow(dat2) # 53

# Drop unnecessary columns
dat2 <- dat2[,c(1,6,7)]

# Change the rownames
rownames(dat2) <- seq(1:53)

# What is the index of the gene the reads were simulated from?
which(dat2$Geneid %in% "ENSG00000151729.11") # 12

# Add a column with the number of simulated reads
simulated <- rep(0,53) 
simulated[12] <- 1000
dat2$simulated_reads <- simulated 

# Write to file 
write.csv(dat2, "output/ensembl_SLC25_featureCounts.csv")
