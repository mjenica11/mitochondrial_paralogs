library(data.table)

# Read in the salmon quant.sf file: contains the number of reads mapped
# to each transcript in the reference transcriptome
dat <- read.csv("SLC25A4_1000_mapped_reads/quant.sf", sep="\t")

# Rename the "Name" column so it can be merged with another dataframe
colnames(dat)[1] <- "RefSeq_ID"
head(dat)

# Subset the SLC25 transcripts from the quant.sf file
# Chose the transcripts that are present in both RefSeq and Ensembl as the
# representative transcript for each SLC25 paralog
# NM_001317231.2 (SLC25A50) is a duplicate of NM_014342.4
RefSeq <-c("NM_032315.3", "NM_207348.3", "NM_013386.5", "NM_014655.4",
		   "NM_003705.5", "NM_017875.4", "NM_000387.6", "NM_173471.4",
		   "NM_001104647.3", "NM_031291.4", "NM_021833.5", "NM_001151.4",
		   "NM_138773.4", "NM_001349336.2", "NM_031947.4", "NM_001271641.2",
		   "NM_004277.5", "NM_018843.4", "NM_014251.3", "NM_016612.4",
		   "NM_030780.5", "NM_033412.4", "NM_001330988.2", "NM_152707.4",
		   "NM_031212.4", "NM_001317231.2", "NM_182556.4", "NM_003355.3",
		   "NM_003356.4", "NM_001191061.2", "NM_002635.4", "NM_014252.4",
		   "NM_001010875.4", "NM_030631.4", "NM_001039355.3", "NM_207117.4",
		   "NM_003562.5", "NM_001320870.2", "NM_001143780.3", "NM_001126121.2",
		   "NM_012140.5", "NM_001034172.4", "NM_173637.4", "NM_024103.3",
		   "NM_001321544.2", "NM_031481.3", "NM_005984.5", "NM_006358.4",
		   "NM_001636.4", "NM_001012755.5", "NM_145305.3", "NM_001152.5",
		   "NM_001282195.2")

# Dataframe of just SLC25 read counts
dat2 <- subset(dat, dat$Name %in% RefSeq)
nrow(dat2) # 53

# Append column with the HGNC names of the SLC25 transcripts for readability
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

# Create a dataframe with the HGNC names and RefSeq IDs
# Add a column with the HGNC gene names
temp <- data.frame(HGNC_names = SLC25_genes, RefSeq_ID = RefSeq) 

# Bind the name dataframe to the quant.sf dataframe and merge by "Name"
# to make sure the rows bind in the same order
# Manually checked to make sure the RefSeq IDs matched the HGNC IDs
dat2 <- merge(dat, temp, by = "RefSeq_ID")
head(dat2) 

# Move the HGNC_names colum to the begining
dat2 <- dat2[,c(6,1,2,3,5,4)]
head(dat2)

# The number of reads simulated per transcript
# Instantiate vector of 0s equal to the number of SLC25 family members present
expected <- rep(0, nrow(dat2))

# Insert 200 at the index for the transcript that the reads were  simulated from (SLC25A4) 
which(dat2$RefSeq_ID %in% "NM_001151.4") #9
expected[9] <- 1000
expected

# Add a column with the number of simulated transcripts
dat2$simulated_fragments <- expected
dat2[9:20,]

# Write to file
write.csv(dat2, "output/1000_reads_SLC25_transcripts.csv")
