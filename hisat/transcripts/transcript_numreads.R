# Extract the SLC25 transcripts from the trans_number_reads.count file

# Read in the file counting the number of reads mapped to each transcript
dat <- read.csv("trans_number_reads.count", sep="\t", header=FALSE)

# Names of the SLC25 transcripts present 
SLC25_genes  <- c("SLC25A33", "SLC25A34", "SLC25A24", "SLC25A44", "SLC25A12",
				  "SLC25A38", "SLC25A20", "SLC25A26", "SLC25A36", "SLC25A31",
				  "SLC25A7", "SLC25A4", "SLC25A46", "SLC25A48", "SLC25A2",
				  "SLC25A49", "SLC25A27", "SLC25A40", "SLC25A13", "SLC25A37",
				  "SLC25A32", "SLC25A51", "SLC25A25", "SLC25A16", "SLC25A28",
				  "SLC25A50", "SLC25A45", "SLC25A8", "SLC25A9", "SLC25A22",
				  "SLC25A3", "SLC25A15", "SLC25A30", "SLC25A21", "SLC25A29",
				  "SLC25A47", "SLC25A11", "SLC25A35", "SLC25A39", "SLC25A19",
				  "SLC25A10", "SLC25A52", "SLC25A41", "SLC25A23", "SLC25A42",
				  "SLC25A18", "SLC25A1", "SLC25A17", "SLC25A6", "SLC25A53",
				  "SLC25A43", "SLC25A5", "SLC25A14")

# Indices for the SLC25 family transcripts in RefSeq annotation
# The HTSeq coutn output asses "rna-" to each of the RefSeq IDs 
RefSeq <-c("rna-NM_032315.3", "rna-NM_207348.3", "rna-NM_013386.5", "rna-NM_014655.4",
		   "rna-NM_003705.5", "rna-NM_017875.4", "rna-NM_000387.6", "rna-NM_173471.4",
		   "rna-NM_001104647.3", "rna-NM_031291.4", "rna-NM_021833.5", "rna-NM_001151.4",
		   "rna-NM_138773.4", "rna-NM_001349336.2", "rna-NM_031947.4", "rna-NM_001271641.2",
		   "rna-NM_004277.5", "rna-NM_018843.4", "rna-NM_014251.3", "rna-NM_016612.4",
		   "rna-NM_030780.5", "rna-NM_033412.4", "rna-NM_001330988.2", "rna-NM_152707.4",
		   "rna-NM_031212.4", "rna-NM_001317231.2", "rna-NM_182556.4", "rna-NM_003355.3",
		   "rna-NM_003356.4", "rna-NM_001191061.2", "rna-NM_002635.4", "rna-NM_014252.4",
		   "rna-NM_001010875.4", "rna-NM_030631.4", "rna-NM_001039355.3", "rna-NM_207117.4",
		   "rna-NM_003562.5", "rna-NM_001320870.2", "rna-NM_001143780.3", "rna-NM_001126121.2",
		   "rna-NM_012140.5", "rna-NM_001034172.4", "rna-NM_173637.4", "rna-NM_024103.3",
		   "rna-NM_001321544.2", "rna-NM_031481.3", "rna-NM_005984.5", "rna-NM_006358.4",
		   "rna-NM_001636.4", "rna-NM_001012755.5", "rna-NM_145305.3", "rna-NM_001152.5",
		   "rna-NM_001282195.2")

# Subset the rows where the SLC25 transcripts are present
dat2 <- dat[is.element(dat$V1, RefSeq),]

# Add a column with the SLC25 names
dat3 <- cbind(data.frame(gene=SLC25_genes), dat2)

# Rename last two columns
colnames(dat3)[2:3] <- c("transcript_ID", "NumReads")

# Write to file
write.csv(dat3, "SLC25_transcript_num_reads.csv")
