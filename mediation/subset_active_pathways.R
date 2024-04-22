# Mediation analysis of SLC carriers
library(medflex)
library(data.table)
library(dplyr)
library(stringr)
library(tidyr)
library(stringr)
library(purrr)

# Read in the heart and liver pathologist activiry scores
activity <- fread("/scratch/mjpete11/mitochondrial_paralogs/linear_models/data/data/Merge_Heart_Liver_pathway_activity.txt", header=TRUE, sep="\t") # float
dim(activity) # 1324 297 
activity[1:5,1:5]

# Make sure the object is the right class
class(activity) # data.table data.frame
activity <- as.data.frame(activity)
class(activity) #  data.frame

# Rename the first column to remove spaces
colnames(activity)[1] <- 'Pathway_Name'

# Store the pathway names to its own object
pathway_names <- as.data.frame(activity$'Pathway_Name')
class(pathway_names) # data.frame 

# Make sure the columns are the correct class
head(apply(activity, 2, class)) # character
act  <- activity[,2:ncol(activity)] %>% mutate_if(is.character, is.numeric)
head(apply(act, 2, class)) # character
head(act)

# Re-append the pathway names column
act1<- cbind(act, pathway_names)
head(act1)

# Re-name the pathway names column
dim(act1) # 1324 297
colnames(act1)[297] <- c('pathway_name')
tail(act1)

# Move the last column to the front
act2 <- act1 %>% dplyr::select(pathway_name, everything())
act2[1:5,1:5]
dim(act2)==dim(activity) # TRUE TRUE
dim(act2) # 1324 297

# Subset pathway names with at least 50% variability in 90% of the samples
df1 <- act2 %>% 
    as_tibble() %>% 
    slice(1:ncol(act2)) %>% 
    filter(any(c_across(starts_with("GTEX")) >= 0.5))

head(df1)
df1[1:5,1:5]
dim(df1) # 297 297
class(df1)

# Don't keep rows with NA values across
library(purrr)
df2 <- df1 %>% filter(if_any(everything(), purrr::negate(is.na))) 

df2 <- df1[complete.cases(df1),]
df2[1:15,1:5]

# Write to file
write.table(df2, 
            file="/scratch/mjpete11/mitochondrial_paralogs/linear_models/data/data/filtered_activity_scores_50percent.csv", 
            sep=",", 
            row.names=FALSE)

################################### Exploratory analysis #########################################################
# Function to filter the dataframe for genes with a certain probability scores in a certain proportion of samples
filter_data <- function(DATA, PROPORTION, PROBABILITY){

    # calculate the proportion of values above a certain probability in each row
    prop_above_threshold <- rowMeans(DATA > PROPORTION)

    # filter the rows where the proportion is >= the probability value in a certain number of samples
    res <- DATA[prop_above_threshold >= PROBABILITY, ]
    return(res)
}

# Strict threshold: How many of the pathways have a probability of being active >90% in at least half the samples
filtered_data <- filter_data(DATA=activity, PROPORTION=0.5, PROBABILITY=0.9)
head(filtered_data)
nrow(filtered_data) # when PROP=0.5, PROB=0.3, nrow=71

# Strict threshold: How many of the pathways have a probability of being active >50% in at least half the samples
filtered_data <- filter_data(DATA=activity, PROPORTION=0.5, PROBABILITY=0.5)
head(filtered_data)
nrow(filtered_data) # when PROP=0.5, PROB=0.3, nrow=139

# Weak threshold: How many of the pathways have a probability of being active >30% in at least half the samples
filtered_data <- filter_data(DATA=activity, PROPORTION=0.5, PROBABILITY=0.3)
head(filtered_data)
nrow(filtered_data) # when PROP=0.5, PROB=0.3, nrow=218

# Try again but set the proportion to 0.9
# Strict threshold: How many of the pathways have a probability of being active >90% in at least half the samples
filtered_data <- filter_data(DATA=activity, PROPORTION=0.9, PROBABILITY=0.9)
head(filtered_data)
nrow(filtered_data) # nrow=6

# Strict threshold: How many of the pathways have a probability of being active >50% in at least half the samples
filtered_data <- filter_data(DATA=activity, PROPORTION=0.9, PROBABILITY=0.5)
head(filtered_data)
nrow(filtered_data) # nrow=9

# Weak threshold: How many of the pathways have a probability of being active >30% in at least half the samples
filtered_data <- filter_data(DATA=activity, PROPORTION=0.9, PROBABILITY=0.3)
head(filtered_data)
nrow(filtered_data) # nrow= 20

# Try a different filtering approach so I can keep the first column
df1 <- activity %>%
    filter(!is.na(.[[1]])) %>% # Keep rows with a value in the first column
    select(-1) %>% # Select all cols except the first
    filter(rowMeans(. >= 0.1) >= 0.5) # Filter on remaining cols
head(df1)
class(df1)
################################### Exploratory analysis #########################################################

# Read in the qnorm + voom normalized heart and liver counts
counts <- fread("/scratch/mjpete11/mitochondrial_paralogs/linear_models/data/data/batch_voom_qnorm_matrix.csv", sep=",")
dim(counts) # 51259 298
counts[1:5,1:5]

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
genes

# Subset the SLC25 carrier and the metabolic biomarkers
class(counts) # data.table data.frame
counts <- as.data.frame(counts)
class(counts) # data.frame
sub_df <- counts[counts$"hugo_names" %in% genes, ]
sub_df

# Which genes MCF genes and/or biomarkers are missing from the filtered/normalized GTEx data
# PFKX, PGAM2, CACT, ID3HA are missing 
setdiff(genes, sub_df$'hugo_names') 

# Write to file
#write.table(sub_df, "/scratch/mjpete11/linear_models/data/biomarkers_sub_df_preprocess.csv", sep=",")

# Read in sample attributes files
file2 <- read.csv("/scratch/mjpete11/mitochondrial_paralogs/linear_models/data/data/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt", sep="\t")

# Subset heart left ventricle and liver samples into two separate dfs
heart <- file2[file2$SMTSD %in% "Heart - Left Ventricle", ]
liver <- file2[file2$SMTSD %in% "Liver", ]

# Subset the SLC25 gene count df by the sample IDs that match the IDs in file3 df
heart2 <- sub_df %>% select(contains(heart$SAMPID))
liver2 <- sub_df %>% select(contains(liver$SAMPID))

# Append the gene ensemble ID and common name
heart3 <- cbind('gene'=sub_df$'hugo_names', heart2)
liver3 <- cbind('gene'=sub_df$'hugo_names', liver2)

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
organs <- organs[!duplicated(organs), ]
head(organs) # subjid gene organ sampid counts smrin smtsisch

# Add the metabolic pathways that were significantly differentially expressed after MaREA 
# to the dataframe
# The number of rows between the reaction_activity_scores and the organs df are not equal
# Add an index to each df so there is a shared column to bind on
organs$index <- seq(from=1, to=nrow(organs), by=1)  
head(organs);tail(organs)

reaction_activity_scores$index <- seq(from=1, to=nrow(reaction_activity_scores), by=1)  
head(reaction_activity_scores);tail(reaction_activity_scores)

# Use merge() to horizontal bind; the shorter df (reaction_activity_scores) will be recycled
combined_data <- merge(organs, reaction_activity_scores, by="index")
head(combined_data);head(combined_data)

# Rename the last column so it is a continuous string
colnames(combined_data)[11] <- "log2_fc"
class(combined_data$log2_fc) # numeric
class(combined_data$SMRIN) # numeric
class(combined_data$integer) # numeric

# Write organs df to file
#write.table(organs, "/scratch/mjpete11/linear_models/data/organs_biomarkers_combat_seq3.csv", sep=",")

# Read organs df back in
#organs <- read.csv("/scratch/mjpete11/linear_models/data/organs_biomarkers_combat_seq3.csv", sep=",")

# How many paired samples (heart and liver from the sample individual) are there?
length(unique(organs$SUBJID)) # 148

# Expand the data
# Use only the SLCs as the predictor variable
#Yimp <- glm(organ ~ gene + SMRIN + SMTSISCH, family="binomial", data=organs, weights=NULL, unlist(subset(organs, grepl("SLC25A",gene), gene))) 

############################# TEST ############################################
# Test if the effect of MCF expression on metabolic pathway state is mediated by metabolic biomarker expression
# Do basic mediation instead of using medflex approch to start
# Split the combined dataset into separate dfs-one for the MCFs and one for the metabolic biomarkers
MCF_df <- combined_data[combined_data$gene %in% SLC,]
metabolic_biomarker_df <- combined_data[combined_data$gene %in% biomarkers,]
# Recombine into one df but add them horizontally instead of vertically
horizontal_df <- merge(x = data.frame(MCF_df, row.names=NULL), 
                       y = data.frame(metabolic_biomarker_df, row.names=NULL), 
                       by = NULL, all.x = TRUE, all.y = TRUE)

# It only kept the heart sample for some reason :(
head(horizontal_df)
tail(horizontal_df)
# Every pathway ~ every gene + rna integrity number + total ischemic time
# The columns that end in .x are the MCF data and the columns that end in .y are the metabolic biomarker data
# Step 1: Estimate the total effect
fit <- lm()

# Step 2: Path A (X on M)

# Step 3:Path B (M on Y, controlling for X)

# Step 4: Reversed Path C (Y on X, controlling for M)

# Summary

################# Make a dataframe out of p-values from the results of a lm() model ##############
# Examine the output
names(Yimp)
Yimp$coefficients
Yimp$residuals
Yimp$effects
summary(Yimp)
names(summary(Yimp))
summary(Yimp)$r.squared # 0.529
summary(Yimp)$adj.r.squared # 0.137
head(summary(Yimp)$coefficients)
# p-values for each independent variable
length(summary(Yimp)$coefficients[,4]) # 120, vector
head(summary(Yimp)$coefficients[,4])
head(summary(Yimp)$coefficients[,4])
# Are any p-values significant?
any(summary(Yimp)$coefficients[,4])<0.05 # FALSE

# Reshape the p-values into a readable table
# Reshape vector to a matrix (single column)
class(summary(Yimp)$coefficients[,4])
pval_df <- as.data.frame(matrix(summary(Yimp)$coefficients[,4], nrow=20, ncol=6), byrow=TRUE)
dim(pval_df)
colnames(pval_df)
head(pval_df)

# Add the gene names as columns to the dataframe
names(summary(Yimp)$coefficients[,4])
length(names(summary(Yimp)$coefficients[,4]))

# Split the gene names vector into separate vectors and then append them as columns
gene_name_vec <- names(summary(Yimp)$coefficients[,4])
gene_name_vec <- str_remove_all(gene_name_vec, "gene")
length(gene_name_vec)==length(summary(Yimp)$coefficients[,4]) # TRUE
gene_name_vec

# Drop the 'gene' string from the items in the vector
vec1 <- gene_name_vec[1:20]
pval_df$gene_names_col1 <- vec1
vec2 <- gene_name_vec[21:40]
pval_df$gene_names_col2 <- vec2
vec3 <- gene_name_vec[41:60]
pval_df$gene_names_col3 <- vec3
vec4 <- gene_name_vec[61:80]
pval_df$gene_names_col4 <- vec4
vec5 <- gene_name_vec[81:100]
pval_df$gene_names_col5 <- vec5
vec6 <- gene_name_vec[101:120]
pval_df$gene_names_col6 <- vec6
# Move the last column to the first placea
class(pval_df) # data.frame
pval_df <- tibble(pval_df)
class(pval_df) # tbl_df tbl data.frame
pval_df <- pval_df %>% relocate(gene_names_col1)
pval_df <- pval_df %>% relocate(gene_names_col2, .after=V1)
pval_df <- pval_df %>% relocate(gene_names_col3, .after=V2)
pval_df <- pval_df %>% relocate(gene_names_col4, .after=V3)
pval_df <- pval_df %>% relocate(gene_names_col5, .after=V4)
pval_df <- pval_df %>% relocate(gene_names_col6, .after=V5)
pval_df <- pval_df %>% relocate(rownames(pval_df), .after=gene_names_col1)
head(pval_df)

# Write the resulting p-values to a table
write.table(pval_df, "/scratch/mjpete11/mitochondrial_paralogs/linear_models/data/data/test_mediation_pvals.csv", sep=",", row.names=FALSE)

############################# TEST ############################################
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

# Test again
dat1 <- subset(organs, gene %in% c("SLC25A5", "MT-ND1"))
dat1$gene <- as.factor(dat1$gene)
dat1$organ <- as.factor(dat1$organ)
dat2 <- dat1[!duplicated(dat1[c(2,3,4)]),]
dat3 <- spread(dat2, key=gene, value=counts)
# Add a column with the pathway state
dat3$pathway_state <- rep("citrate_cycle", nrow(dat3))
dat4 <- dat3[complete.cases(dat3),]
head(dat4)
########################### Mediation analysis ########################### 

# Function to make list of dataframes for mediation
make_dfs <- function(GENE, BIOMARKER){
	dat <- subset(organs, gene %in% c(GENE, BIOMARKER)) # Subset the genes of interest
	dat$gene <- as.factor(dat$gene) # Convert to factor
	dat$organ <- as.factor(dat$organ) # Convert to factor 
	#dat <- dat[!duplicated(dat$SAMPID),] # Drop duplicate rows
	dat <- dat[!duplicated(dat[c(2,3,4)]),] # Drop duplicate rows

	# Reshape test dataframe; make SLC and biomarker genes into seperate columns
	dat2 <- spread(dat, key=gene, value=counts)
	dat3 <- dat2[complete.cases(dat2),] # Drop NAs
	dat3 # Print so I can check the data frame
	return(dat3)
}
#dats <- Map(make_dfs, GENE=SLC[1])
result1 <- make_dfs(GENE=SLC[1], BIOMARKER=biomarkers[7])
head(result1)

############## By hand since sometime I get singularity errors ###############
#rm(list=c("dat1","Yimp0","impData0","Yfit0","res"))
#which(biomarkers=="MT-ND1")
#which(SLC=="SLC25A5")
#dat1 <- make_dfs(GENE=SLC[6], BIOMARKER=biomarkers[70])
# Rename MT-DN1 column because the glm function can't detect dashes
#colnames(dat4)[6] <- "MTND1"
#head(dat4)

############################## Mediation ######################################
rm(list=c("dat1","Yimp0","impData0","Yfit0","res"))
# Expand data: Step 1: fit glm
Yimp0 <- glm(pathway_state ~ SLC25A1 + PFKM + SMRIN + SMTSISCH, family="binomial", data=result1) 
attributes(Yimp0)
head(Yimp0[["effects"]])
# Expand data: Step 2: impute data 
impData0 <- neImpute(Yimp0)
head(impData0)
# Fit the natureal effect model
Yfit0 <- neModel(organ ~ SLC25A10 + SLC25A11, expData=impData0, se="robust", family="binomial") 
# Model fit summary
summary(neEffdecomp(Yfit0)) # Get the natural direct and indirect effect
res <- summary(neEffdecomp(Yfit0)) # Get the total effect
# Write plot to file
pdf("SLC25A1_PFKM.pdf")
plot(neEffdecomp(Yfit0), xlab="Effest size estimate", main="SLC25A1 (A) and PFKM (M)")
dev.off()
# Write summary to file
sink("SLC25A1_PFKM.txt")
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

