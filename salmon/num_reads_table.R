library(data.table)

# Read in the quant.sf files
files <- list.files("individual_simulations", pattern="*.sf", recursive=TRUE, full.names=TRUE)
res <- lapply(files, function(i){
					  read.csv(i, header=TRUE, sep='\t')
})

# Rename the "Name" column so it can be merged with another dataframe
colnames <- c("RefSeq_ID", "Length", "EffectiveLength", "TPM", "NumReads")
for (i in seq_along(res)){
		colnames(res[[i]]) <- colnames 
}

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

# HGNC names of SLC25 paralogs 
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

# List of genes that are not paralogs as negative controls
not_paralogs <- c("NR_001564.2","NM_001385731.1", "NM_014421.3", "NM_030636.3",
				  "NM_033082.4", "NM_019598.3", "NM_001271885.2", "NM_001077503.2", 
				  "NM_001281.3", "NM_003489.4")

not_paralog_names <- c("XIST", "CSN2", "DKK2", "EEPD1", "SARNP", "YLMP1", "AAGAB",
					   "TAC4", "TBCB", "NRIP1")

# Combine list of SLC25 genes and non-paralog controls
all_RefSeq <- c(RefSeq, not_paralogs)
all_genes <- c(SLC25_genes, not_paralog_names)

# Dataframe of genes that are not members of the SLC25 family
'%!in%' <- Negate('%in%')
not <- lapply(res, function(x) subset(x, x$RefSeq_ID %!in% RefSeq))

# Check if any rows have counts > 0
#tmp <- lapply(not, function(x) print(x[x$NumReads > 0])) # None

# Dataframe of SLC25 read counts and non-paralogous controls
dat <- lapply(res, function(x) subset(x, x$RefSeq_ID %in% all_RefSeq))
lapply(dat, nrow)

# Create a dataframe with the HGNC names and RefSeq IDs
# Add a column with the HGNC gene names
temp <- data.frame(HGNC_names = all_genes, RefSeq_ID = all_RefSeq) 

# Bind the name dataframe to the quant.sf dataframe and merge by "RefSeq_ID"
# to make sure the rows bind in the same order
# Manually checked to make sure the RefSeq IDs matched the HGNC IDs
dat2 <- lapply(dat, function(x) merge(x, temp, by="RefSeq_ID"))
head(dat2[[1]]);head(dat2[[53]])

# Move the HGNC_names column to the begining
dat2 <- lapply(dat2, function(x) x[,c(6,1,2,3,5,4)])

# The number of reads simulated per transcript
expected <- rep(1000, nrow(dat2[[1]]))
expected
dat3 <- lapply(dat2, function(x) cbind(x, simulated_fragments = expected))

# Combine each df into one
genes <- c("SLC25A1","SLC25A10","SLC25A11","SLC25A12","SLC25A13","SLC25A14",
		  "SLC25A15","SLC25A16","SLC25A17","SLC25A18","SLC25A19","SLC25A2",
		  "SLC25A20","SLC25A21","SLC25A22","SLC25A23","SLC25A24","SLC25A25",
		  "SLC25A26","SLC25A27","SLC25A28","SLC25A29","SLC25A3","SLC25A30",
		  "SLC25A31","SLC25A32","SLC25A33","SLC25A34","SLC25A35","SLC25A36",
		  "SLC25A37","SLC25A38","SLC25A39","SLC25A4","SLC25A40","SLC25A41",
		  "SLC25A42","SLC25A43","SLC25A44","SLC25A45","SLC25A46","SLC25A47",
		  "SLC25A48","SLC25A49","SLC25A5","SLC25A50","SLC25A51","SLC25A52",
		  "SLC25A53", "SLC25A6", "SLC25A7", "SLC25A8", "SLC25A9")

# Function to slice relevant row from each df in list of dataframes
slice_dataframes <- function(vec, name){
		dx <- dat3[[vec]][which(dat3[[vec]]$HGNC_names %in% name),]
		return(dx)
}

# Apply function to list of dfs
df_lst <- Map(slice_dataframes, vec = seq(1:53), name = genes)
df_lst

# Subset non-paralogs from each df in list
slice2 <- function(vec) {
	dy <- dat3[[vec]][which(dat3[[vec]]$HGNC_names %in% not_paralog_names),]
	return(dy)
}
df_lst2 <- Map(slice2, vec = seq(1:53))

# Do any of the non-paralogous genes have counts > 0?
res2 <- lapply(df_lst2, function(x) x[which(x$NumReads > 0),]) # No

# Merge dataframes
cols <- c("HGNC_names", "RefSeq_ID", "Length", "EffectiveLength", "NumReads", "TPM", "simulated_fragments")
dat4 <- Reduce(function(...) merge(..., all = TRUE), df_lst)
head(dat4)
tail(dat4)

# Append one of the non-paralogous dataframes to the list of paralogous genes
tmp <- rbind(dat4, res2[[1]])

# Add a column with the true positive rate
dat4$true_positive <- dat4$NumReads/dat4$simulated_fragments * 100
min(dat4$true_positive) 
max(dat4$true_positive) 

# Write to file
write.csv(dat4, "output/1000_reads_individual_transcripts.csv")
