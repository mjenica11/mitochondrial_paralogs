# Mediation analysis of SLC carriers
library(medflex)
library(data.table)
library(dplyr)
library(stringr)
library(tidyr)

# Read in preprocessCore quantile normalized counts
counts <- fread("/scratch/mjpete11/linear_models/data/voom_quantile_normalized_counts.csv", sep=",") # float

# Drop the index column
counts$V1 <- NULL

# List of SLC25 paralogs
SLC <- c("SLC25A1", "SLC25A2", "SLC25A3", "SLC25A4", "SLC25A5", "SLC25A6",
		 "UCP1", "UCP2", "UCP3", "SLC25A10", "SLC25A11", "SLC25A12", 
		 "SLC25A13", "SLC25A14", "SLC25A15", "SLC25A16", "SLC25A17", 
		 "SLC25A18", "SLC25A19", "SLC25A20", "SLC25A21", "SLC25A22", 
		 "SLC25A23", "SLC25A24", "SLC25A25", "SLC25A26", "SLC25A27", 
		 "SLC25A28", "SLC25A29", "SLC25A30", "SLC25A31", "SLC25A32", 
		 "SLC25A33", "SLC25A34", "SLC25A35", "SLC25A36", "SLC25A37", 
		 "SLC25A38", "SLC25A39", "SLC25A40", "SLC25A41", "SLC25A42", 
		 "SLC25A43", "SLC25A44", "SLC25A45", "SLC25A46", "SLC25A47",
		 "SLC25A48", "MTCH1", "MTCH2", "SLC25A51", "SLC25A52", "SLC25A53")

# List of metabolic biomarkers 
biomarkers <- c("SLC16A1", "SLC16A2", "HK1", "HK2", "HK3", "GPI", "PFKM",
				"PFKX", "PFKP", "ALDOA", "ALDOB", "TPI1", "GAPDH", "PGK1",
				"PGK2", "PGAM1", "PGAM5", "PGAM2", "PGAM4", "ENO1",
				"PKLR", "LDHA", "LDHB", "LDHC", "LDHD", "CACT", "CPT1A",
				"CPT1B", "CPT2", "ACADVL", "ACADL", "ACADM", "ACADS",
				"ACADSB", "ACAD11", "ACAD8", "ACAD10", "ACAD9", "ECHS1",
				"ECH1", "ECI1", "ECI2", "ECHS1", "HADHA", "PDHA1", "PDHA2",
				"PDHB", "PDP1", "PDHX", "CS", "ACO1", "ACO2", "IDH3G",
				"IDH1", "IDH2", "ID3HA", "IDH3B", "OGDH", "SUCLG2",
				"SUCLG1", "SUCLA2", "SDHB", "SDHD", "SDHA", "SDHC", "SDHAF2",
				"FH", "MDH1", "MDH2", "MT-ND1")

genes <- c(SLC, biomarkers)

# Subset the SLC25 carrier and the metabolic biomarkers
counts <- as.data.frame(counts)
sub_df <- counts[counts$"Description" %in% genes, ]

# PFKX, PGAM2, CACT, ID3HA are missing 
setdiff(genes, sub_df$'Description') 

# Write to file
write.table(sub_df, "/scratch/mjpete11/linear_models/data/biomarkers_sub_df_preprocess.csv", sep=",")

# Read in sample attributes files
file2 <- read.csv("/scratch/mjpete11/linear_models/data/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt", sep="\t")

# Subset heart left ventricle and liver samples into two separate dfs
heart <- file2[file2$SMTSD %in% "Heart - Left Ventricle", ]
liver <- file2[file2$SMTSD %in% "Liver", ]

# Subset the SLC25 gene count df by the sample IDs that match the IDs in file3 df
heart2 <- sub_df %>% select(contains(heart$SAMPID))
liver2 <- sub_df %>% select(contains(liver$SAMPID))

# Append the gene ensemble ID and common name
heart3 <- cbind('gene'=sub_df$'Description', heart2)
liver3 <- cbind('gene'=sub_df$'Description', liver2)

# Add a columns with the organ in each df
heart3$organ <- 'heart' 
liver3$organ <- 'liver' 

# Reshape dataframe so it can be converted to a design matrix object
heart4 <- reshape2::melt(data = heart3,
			   id.vars = c("gene", "organ"),
			   measure.vars = colnames(heart3)[3:ncol(heart3)-1],
			   variable.name = "samples",
			   value.name = "counts")

liver4 <- reshape2::melt(data = liver3,
			   id.vars = c("gene", "organ"),
			   measure.vars = colnames(liver3)[3:ncol(liver3)-1],
			   variable.name = "samples",
			   value.name = "counts")

colnames(heart4)[3] <- "SAMPID" 
colnames(liver4)[3] <- "SAMPID" 

# Add column with the expression batch ID (SMGEBTCH) and the type of genotype
# or expression batch (SMGEBTCHT)
heart4$SMGEBTCH <- file2$SMGEBTCH[match(heart4$SAMPID, file2$SAMPID)]
liver4$SMGEBTCH <- file2$SMGEBTCH[match(liver4$SAMPID, file2$SAMPID)]

heart4$SMGEBTCHT <- file2$SMGEBTCHT[match(heart4$SAMPID, file2$SAMPID)]
liver4$SMGEBTCHT <- file2$SMGEBTCHT[match(liver4$SAMPID, file2$SAMPID)]

heart4$SMRIN <- file2$SMRIN[match(heart4$SAMPID, file2$SAMPID)]
liver4$SMRIN <- file2$SMRIN[match(liver4$SAMPID, file2$SAMPID)]

heart4$SMTSISCH <- file2$SMTSISCH[match(heart4$SAMPID, file2$SAMPID)]
liver4$SMTSISCH <- file2$SMTSISCH[match(liver4$SAMPID, file2$SAMPID)]

# Combine into one df
organs <- rbind(heart4, liver4)

# Add sample ID col to meta
heart4$SUBJID <- str_sub(heart4$SAMPID, start=1L, end=10L)
liver4$SUBJID <- str_sub(liver4$SAMPID, start=1L, end=10L)

# Subset to only heart and liver samples that are from the same person
heart4$SMGEBTCH <- NULL
#heart4$SMRIN <- NULL
#heart4$SMTSISCH <- NULL
heart4$SMGEBTCHT <- NULL

liver4$SMGEBTCH <- NULL
#liver4$SMRIN <- NULL
#liver4$SMTSISCH <- NULL
liver4$SMGEBTCHT <- NULL

# Dataframe with only samples from people who donated both a heart and liver: organs
# e.g. paired samples
organs_tmp <- merge(heart4, liver4, by="SUBJID") 
length(unique(organs_tmp$SUBJID)) # 148 individuals
tmp1 <- organs_tmp[,c("SUBJID", "gene.x", "organ.x", "SAMPID.x", "counts.x", "SMRIN.x", "SMTSISCH.x")]
tmp2 <- organs_tmp[,c("SUBJID", "gene.y", "organ.y", "SAMPID.y", "counts.y", "SMRIN.y", "SMTSISCH.y")]
colnames(tmp1) <- c("SUBJID", "gene", "organ", "SUBJID", "counts", "SMRIN", "SMTSISCH")
colnames(tmp2) <- c("SUBJID", "gene", "organ", "SUBJID", "counts", "SMRIN", "SMTSISCH")
organs <- rbind(tmp1, tmp2)
colnames(organs)[4] <- "SAMPID"

# Write organs df to file
write.table(organs, "/scratch/mjpete11/linear_models/data/organs_biomarkers_combat_seq3.csv", sep=",")

# Read organs df back in
organs <- read.csv("/scratch/mjpete11/linear_models/data/organs_biomarkers_combat_seq3.csv", sep=",")

# How many paired samples (heart and liver from the sample individual) are there?
length(unique(organs$SUBJID)) # 148

########################### Mediation analysis ########################### 
# Expand the data
# Use only the SLCs as the predictor variable
#Yimp <- glm(organ ~ gene + SMRIN + SMTSISCH, family="binomial", data=organs, weights=NULL, unlist(subset(organs, grepl("SLC25A",gene), gene))) 

############################# TEST ############################################
#Yimp <- glm(organ ~ gene + SMRIN + SMTSISCH, family="binomial", data=organs, weights=NULL, unlist(subset(organs, grepl("SLC25A",gene), gene))) 

# Subset df for testing and replace some of the heart vars with liver 
# so there are 2 factor levels for the organs var
#tmp_organs <- organs[100:140,]
#tmp_organs$organ[1:10] <- "liver"
#tmp_organs$SUBJID <- NULL # Delete unused columns or else glm() complains
#tmp_organs$SAMPID <- NULL
tmp_organs <- subset(organs, gene %in% c("SLC25A4", "PFK"))
tmp_organs$gene <- as.factor(tmp_organs$gene)
tmp_organs$organ <- as.factor(tmp_organs$organ) # Convert to factor 
tmp_organs <- tmp_organs[!duplicated(tmp_organs),]
organs <- organs[!duplicated(organs),]

# Reshape test dataframe; make SLC and biomarker genes into seperate columns
tmp3 <- spread(tmp_organs, key=gene, value=counts)
head(tmp3)

# Function to make list of dataframes for mediation
make_dfs <- function(GENE, BIOMARKER){
	dat <- subset(organs, gene %in% c(GENE, BIOMARKER)) # Subset the genes of interest
	dat$gene <- as.factor(dat$gene) # Convert to factor
	dat$organ <- as.factor(dat$organ) # Convert to factor 
	dat <- dat[!duplicated(dat$SAMPID),] # Drop duplicate rows

	# Reshape test dataframe; make SLC and biomarker genes into seperate columns
	dat2 <- spread(dat, key=gene, value=counts)
	dat3 <- dat2[complete.cases(dat2),] # Drop NAs
	dat3 # Print so I can check the data frame
	return(dat3)
}
#dats <- Map(make_dfs, GENE=SLC[1])

dat1 <- subset(organs, gene %in% c("SLC25A6", "MT-ND1"))
dat1$gene <- as.factor(dat1$gene)
dat1$organ <- as.factor(dat1$organ)
dat2 <- dat1[!duplicated(dat1[c(2,3,4)]),]
dat3 <- spread(dat1, key=gene, value=counts)
dat4 <- dat3[complete.cases(dat3),]
head(dat4)

############## By hand since sometime I get singularity errors ###############
rm(list=c("dat1","Yimp0","impData0","Yfit0","res"))
which(biomarkers=="MT-ND1")
which(SLC=="SLC25A6")
dat1 <- make_dfs(GENE=SLC[6], BIOMARKER=biomarkers[70])
# Rename MT-DN1 column because the glm function can't detect dashes
colnames(dat1)[6] <- "MTND1"
head(dat1)
# Expand data: Step 1: fit glm
Yimp0 <- glm(organ ~ MTND1 + SLC25A6 + SMRIN + SMTSISCH, family="binomial", data=dat1) 
attributes(Yimp0)
head(Yimp0[["effects"]])
# Expand data: Step 2: impute data 
impData0 <- neImpute(Yimp0)
head(impData0)
# Fit the natureal effect model
Yfit0 <- neModel(organ ~ MTND10 + MTND11, expData=impData0, se="robust", family="binomial") 
# Model fit summary
summary(neEffdecomp(Yfit0)) # Get the natural direct and indirect effect
res <- summary(neEffdecomp(Yfit0)) # Get the total effect
# Write plot to file
pdf("MTND1_SLC25A6.pdf")
plot(neEffdecomp(Yfit0), xlab="Effest size estimate", main="MTND1 (A) and SLC25A6 (M)")
dev.off()
# Write summary to file
sink("MTND1_SLC25A6.txt")
print(res)
sink()
############## By hand since sometime I get singularity errors ###############



# Odds ratio of the natural indirect effect
# Example: Changing the biomarker ACO1 value that would have been observed at
# the low level of SLC25A11 (M(0)) to the value that would have been observed
# at high level of SLC25All (M(1)) while keeping SLC high (A=1), increases
# the odds of the observed outcome variable (organ) by...
exp(-9.259e-04) # 1

############################# Mediation Function  ############################################
# Mediaton function
mediation <- function(GENE, NUM){
	dat <- subset(organs, gene %in% c(GENE, "HK1")) # Subset the genes of interest
	dat$gene <- as.factor(dat$gene) # Convert to factor
	dat$organ <- as.factor(dat$organ) # Convert to factor 
	dat <- dat[!duplicated(dat),] # Drop duplicate rows

	# Reshape test dataframe; make SLC and biomarker genes into seperate columns
	dat2 <- spread(dat, key=gene, value=counts)
	dat3 <- dat2[complete.cases(dat2),] # Drop NAs
	dat3 # Print so I can check the data frame

	# Expand data: Step 1: fit glm
	Yimp <- glm(organ ~ reformulate(GENE) + HK1 + SMRIN + SMTSISCH, family="binomial", data=dat3) 

	# Expand data: Step 2: impute data 
	impData <- neImpute(Yimp)

	# Fit the natureal effect model
	formula0 <- paste(paste0(SLC,"0"),paste0(SLC,"1"), sep=" + ") # Make formulas as strings
	Yfit <- neModel(organ ~ reformulate(paste(formula[NUM])), expData=impData, se="robust", family="binomial") 

	# Model fit summary
	res <- summary(neEffdecomp(Yfit))
	return(res)
}

# Apply mediator function with each SLC as the independent variable
# and hexokinase 1 as the mediator
RES <- Map(mediation, GENE=SLC, NUM=seq(1,53,by=1)) 
RES <- Map(mediation, GENE=SLC[4]) 

###### Some test code
form1 <- paste(GENE0, "HK1", "SMRIN", "SMTSISCH", sep=" + ")
Yimp0 <- glm(organ ~ reformulate(form1), family="binomial", data=tmp3) 
formula0 <- paste(paste0(SLC,"0"),paste0(SLC,"1"), sep=" + ") # Make formulas as strings
Yfit <- neModel(organ ~ reformulate(paste(formula[NUM])), expData=impData, se="robust", family="binomial") 

