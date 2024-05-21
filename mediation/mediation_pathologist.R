# Script to do mediation analysis on the activity scores 

# Library
library(data.table)
library(janitor)
library(tidyverse)
#library(dplyr)
library(stringr)
library(mediation)
library(purrr)
library(broom)

# Read in the pathways and activity scores at 50% activity
activity <- fread("/scratch/mjpete11/mitochondrial_paralogs/linear_models/data/data/filtered_activity_scores_50percent.csv",sep=",")
activity[1:5,1:5]

# Read in the gene counts with the RIN and activity scores
organs <- read.csv("/scratch/mjpete11/mitochondrial_paralogs/linear_models/data/data/organs_biomarkers_combat_seq3.csv", sep=",")
organs[1:5,1:5]

# Change the .s to  _ in the sample names in the activity dataframe
colnames(activity) <- gsub(x=colnames(activity), pattern="\\.", replacement="-")

# Transpose the activity df 
activity <- t(activity)

# Move the first row to the column name
act <- activity %>% row_to_names(row_number=1)
act[1:5,1:5]
tail(act)
class(act) # matrix array
act <- as.data.frame(act)
class(act) # data.frame
lapply(act, class)
dat <- mutate_all(act, function(x) as.numeric(x))
all(lapply(dat, class)=="numeric") # TRUE
head(dat)
tail(dat)
class(dat)
dim(dat) # 296 292
any(as.logical((sapply(dat, is.na))))==TRUE # FALSE 

# Change the rownames to be the first column
dat$SAMPID <- rownames(dat)

# Move the last column to the first
dat <- dat %>% select(SAMPID, everything())
dat[1:5,1:5]

# Drop the rownames column
rownames(dat) <- NULL
dat[1:5,1:5]

# Add the dativity score per sample to the organs df by merging the dataframes
dat1 <- merge(dat, organs, by="SAMPID")
dat1[1:5,1:5]
tail(dat1)
dim(dat1) # 4191656 
dim(dat) #296 292
dim(organs) # 4191656 7

# Drop the duplicate rows
any(duplicated(dat1)) # TRUE
dat2 <- dat1[!duplicated(dat1),]
any(duplicated(dat2)) # FALSE 
dim(dat2) # 35224 298
dat2[1:5,1:5]
tail(dat2)
class(dat2) # data.frame

# Check if a pathway is present in the activity scores list
pathways <- colnames(act)
pathways
res <- str_detect(pathways, "citrate")
which(res, arr.ind=TRUE) # 231
pathways[231] # "citrate cycle (tca cycle)(kegg)"

# Clear the workspace
stuff <- c("act", "activity", "dat", "dat1", "organs", "pathways", "res") 
rm(list=stuff)
ls() # "dat2"  "dat3"  "stuff"

# Check for equality
dim(dat2)==dim(dat3) # TRUE TRUE
# Write the results to file and read them back in so I don't have to re-run all the prevous code
#write.csv(dat2, 
 #         file="/scratch/mjpete11/mitochondrial_paralogs/linear_models/data/data/activity_scores50_and_SLC_and_biomarker_counts.csv",
  #        sep=",",
   #       row.names=FALSE)

######################################################## Mediation analysis ###########################################################
# Read the pathways with 50% activity and the counts/RIN/ischemic values for the SLCs + Biomarkers
dat3 <- read.csv(file="/scratch/mjpete11/mitochondrial_paralogs/linear_models/data/data/activity_scores50_and_SLC_and_biomarker_counts.csv",
                 sep=",")
# Inspect df
dat3[1:5,1:5] # Since there were spaces in the pathway names, write.csv() replaced those with dots

# Split the dataframe with the SLC25 + the biomarker counts + pathway activity scores into 
# two separate dfs
slc25_df <- subset(dat3, grepl('SLC25', x=dat3$gene))
class(slc25_df) # df
dim(slc25_df) # 14503 298
head(slc25_df$gene); tail(slc25_df$gene)
slc25_df[1:5,1:5]

biomarkers_df <- subset(dat3, !grepl('SLC25', x=dat3$gene))
class(biomarkers_df) # df
dim(biomarkers_df) # 20720 298
head(biomarkers_df$gene); tail(biomarkers_df$gene)
biomarkers_df[1:5,1:5]

# Bind dataframes horizontally for the third step in the model (full model)
horizontal_df <- dplyr::bind_rows(biomarkers_df, slc25_df)
dim(horizontal_df)#35224 298 
horizontal_df[1:5,1:5]
tail(horizontal_df)
head(horizontal_df$gene);tail(horizontal_df$gene)

# Split the dataframes by organ
class(horizontal_df$organ) # charater
horizontal_df$organ <- as.factor(horizontal_df$organ)
class(horizontal_df$organ) # factor

split_df <- split(horizontal_df, horizontal_df$organ) # This produces a list of dfs
heart_df <- split_df[[1]]
liver_df <- split_df[[2]]

# Which column is the citrate pathway in the column bind df
head(horizontal_df[,232]) # citratee.cycle..tca.cycle..kegg.

# DOuble check that the order of the colnames are the same
all(colnames(slc25_df)==colnames(biomarkers_df)) # TRUE

# How to subset the count value based on the string match in another column
slc25A1_counts <- subset(horizontal_df$counts, grepl('SLC25A1', horizontal_df$gene, ignore.case=TRUE)) # ..This return the row values, not the cell :(
head(slc25A1_counts)

# Test that I can subset just the biomarkers from the merged dataframe
head(subset(horizontal_df, !grepl('SLC25', x=horizontal_df$gene, ignore.case=TRUE))$gene)
head(subset(horizontal_df, !grepl('SLC25', x=horizontal_df$gene, ignore.case=TRUE))$counts)
class(subset(horizontal_df, !grepl('SLC25', x=horizontal_df$gene, ignore.case=TRUE))$counts) # numeric
tail(colnames(horizontal_df))

# Try a different subset method; YAY! It worked :D
lm(horizontal_df[,232] ~ counts, data=horizontal_df, subset=grepl('SLC25', x=horizontal_df$gene))

############################################### Linear Models ###########################################
# List of SLC25 paralogs
SLC <- c("SLC25A1", "SLC25A2", "SLC25A3", "SLC25A4", "SLC25A5", "SLC25A6","UCP1", "UCP2", "UCP3", "SLC25A10", "SLC25A11", "SLC25A12","SLC25A13", "SLC25A14", "SLC25A15", "SLC25A16", "SLC25A17","SLC25A18", "SLC25A19", "SLC25A20", "SLC25A21", "SLC25A22","SLC25A23", "SLC25A24", "SLC25A25", "SLC25A26", "SLC25A27","SLC25A28", "SLC25A29", "SLC25A30", "SLC25A31", "SLC25A32","SLC25A33", "SLC25A34", "SLC25A35", "SLC25A36", "SLC25A37","SLC25A38", "SLC25A39", "SLC25A40", "SLC25A41", "SLC25A42","SLC25A43", "SLC25A44", "SLC25A45", "SLC25A46", "SLC25A47","SLC25A48", "MTCH1", "MTCH2", "SLC25A51", "SLC25A52", "SLC25A53")

# List of metabolic biomarkers 
biomarkers <- c("SLC16A1", "SLC16A2", "HK1", "HK2", "HK3", "GPI", "PFKM", "PFKP", "ALDOA", "ALDOB", "TPI1", "GAPDH", "PGK1","PGAM1", "PGAM5", "PGAM4", "ENO1","PKLR", "LDHA", "LDHB", "LDHC", "LDHD", "CPT1A","CPT1B", "CPT2", "ACADVL", "ACADL", "ACADM", "ACADS","ACADSB", "ACAD11", "ACAD8", "ACAD10", "ACAD9", "ECHS1","ECH1", "ECI1", "ECI2", "ECHS1", "HADHA", "PDHA1", "PDHA2","PDHB", "PDP1", "PDHX", "CS", "ACO1", "ACO2", "IDH3G","IDH1", "IDH2", "IDH3B", "OGDH", "SUCLG2","SUCLG1", "SUCLA2", "SDHB", "SDHD", "SDHA", "SDHC", "SDHAF2","FH", "MDH1", "MDH2", "MT-ND1")

# Combine list
genes <- c(SLC, biomarkers)
genes
length(genes) # 118

# Function to iterate through the list of mitochondrial carriers and biomarkers
linear_models <- function(GENE){
    res <- lm(heart_df[,232] ~ counts, data=heart_df, subset=grepl(GENE, x=heart_df$gene), na.action=na.exclude)
    return(res)
}
linear_models(GENE=genes)
heart_lms <- Map(linear_models, GENE=genes) 
#heart_summs <- Map(summary, heart_lms)
class(heart_summs[[1]])
heart_summs[[1]]
length(heart_summs) # 118

# Make the list of summary results into a readable table
glance(heart_lms[[1]]) # produces a "tbl_df"     "tbl"        "data.frame"
heart_dfs <- Map(glance, heart_lms)
class(heart_dfs) # list of dfs
head(heart_dfs[[1]])

# Combine list of dfs into one dataframe
heart <- bind_rows(heart_dfs)
head(heart)

# Add a column with the gene name and move it to the front
heart$gene <- names(heart_dfs)
heart <- heart %>% dplyr::select(gene, everything())

# Write to file
write.table(heart, "/scratch/mjpete11/mitochondrial_paralogs/linear_models/data/data/heart_linear_model_summaries.txt", sep="\t", row.names=FALSE)

# Function to iterate through the list of mitochondrial carriers and biomarkers
#It wouldn't accept a dataframe as a function parameter, so I am hard-coding the dataframe
linear_models <- function(GENE){
    res <- lm(liver_df[,232] ~ counts, data=liver_df, subset=grepl(GENE, x=liver_df$gene), na.action=na.exclude)
    return(res)
}
linear_models(GENE=genes)
liver_lms <- Map(linear_models, GENE=genes) 

# Make the list of summary results into a readable table
glance(liver_lms[[1]]) # produces a "tbl_df"     "tbl"        "data.frame"
liver_dfs <- Map(glance, liver_lms)
class(liver_dfs) # list of dfs
head(liver_dfs[[1]])

# Combine list of dfs into one dataframe
liver <- bind_rows(liver_dfs)
head(liver)

# Add a column with the gene name and move it to the front
liver$gene <- names(liver_dfs)
liver <- liver %>% dplyr::select(gene, everything())

# Write to file
write.table(liver, "/scratch/mjpete11/mitochondrial_paralogs/linear_models/data/data/liver_linear_model_summaries.txt", sep="\t", row.names=FALSE)

# Tried using Map() but I am getting a weird error: Error in lm.fit(x, y, offset = offset, singular.ok = singular.ok, ...) :                                                                                    
#  0 (non-NA) cases --> implies that one of the dependent variables has the same values as the independent variables                       

# Try to figure out which gene is causing the error message
#Map(linear_models, GENE=genes[[83]]) 
# item 61 : "PKKX"; I think this was a typo and it is not a gene
# item 70 : "PGAM2"
# item 70 : "CACT"
# item 106 : "ID3HA"
# item  : "PGAM2"

# Function to iterate through the list of mitochondrial carriers and biomarkers; add technial covariates
linear_models_technical_covariates <- function(GENE){
    res <- lm(liver_df[,232] ~ counts + SMRIN + SMTSISCH, data=liver_df, subset=grepl(GENE, x=liver_df$gene), na.action=na.exclude)
    return(res)
}
liver_cov_lms <- Map(linear_models_technical_covariates, GENE=genes) 

# Make the list of summary results into a readable table
glance(liver_cov_lms[[1]]) # produces a "tbl_df"     "tbl"        "data.frame"
liver_cov_dfs <- Map(glance, liver_cov_lms)
class(liver_cov_dfs) # list of dfs
head(liver_cov_dfs[[1]])

# Combine list of dfs into one dataframe
liver_cov <- bind_rows(liver_cov_dfs)
head(liver_cov)

# Add a column with the gene name and move it to the front
liver_cov$gene <- names(liver_cov_dfs)
liver_cov <- liver_cov %>% dplyr::select(gene, everything())

# Write to file
write.table(liver_cov, "/scratch/mjpete11/mitochondrial_paralogs/linear_models/data/data/liver_technical_covariates_linear_model_summaries.txt", sep="\t", row.names=FALSE)


############################################# Mediation analysis draft #####################################
# Draft approach; skip adding the technical covariates
# Every pathway ~ every gene + rna integrity number + total ischemic time

# Step 1: Estimate the total effect (X on Y); SLC25 genes on the 50% active pathway scores
# The 231st column is the pathway activity scores for the citrate/TCA cycle (kegg) 
model_direct <- lm(slc25_df[,232] ~ slc25_df$gene, data=slc25_df)
summary(model_direct)

# Step 2: Path A (X on M); Estimate the effect of SLC on biomarker
model_mediate <- lm(subset(horizontal_df, !grepl('SLC25', x=horizontal_df$gene, ignore.case=TRUE))$count ~ subset(horizontal_df, grepl('SLC25', x=horizontal_df$gene, ignore.case=TRUE))$count, data=horizontal_df)
model_mediate
names(model_mediate)
model_mediate$effects
names(model_mediate$effects) # I will have to manually figure out how each variable combined together
length(model_mediate$residuals) # 20720 variables
summary(model_mediate)

# Step 3:Path B (M on Y, controlling for X); Estimate the effect of  SLC + Biomarker on pathway
model_indirect <- lm(horizontal_df[,232] ~ counts, data=horizontal_df, subset=grepl('SLC25', x=horizontal_df$gene))
summary(model_indirect)

# Is the mediation effect statistically significant?
results <- mediate(model_mediate, model_indirect, treat='horizontal_df[,232]', mediator='counts', boot=TRUE, sims=500)


############################### Mediation analysis on individual genes ################################################
# Citrate carrier pathway
# First, get the dataframe in the right shape
subdf <- horizontal_df %>% split(grepl('PFKM', x=.$gene))
head(subdf[[2]])

subdf1 <- horizontal_df %>% split(grepl('\\bSLC25A1\\b', x=.$gene)) # Check the subset dfs are the same nrow
head(subdf1[[2]])

nrow(subdf[[2]])==nrow(subdf[[2]]) # TRUE

# Rename the counts column in the first df so I can bind it to the second one
colnames(subdf[[2]])[296] <- "PFKM_counts" 
tail(subdf[[2]])

# Add this column to the subdf1 column
combo_df <- bind_rows(subdf[[2]], subdf1[[2]])
combo_df[1:5,1:5]
tail(combo_df)

# For some reason it put NA in the PFKM counts so I am going to force write it
combo_df$PFKM_counts <- subdf1[[2]]$counts
tail(combo_df)
combo_df$PFKM_counts

# Rename the counts column in the first df so I can bind it to the second one
colnames(combo_df)[299] <- "SLC25A1_counts" 
combo_df$SLC25A1_counts
class(combo_df)

length(combo_df$PFKM_counts)==length(combo_df$SLC25A1_counts) #TRUE

# Step 1: Estimate total effect (x on Y) ; citrate carrier on citrate carrier pathway   
model_direct1 <- lm(combo_df[,232] ~ SLC25A1_counts, data=combo_df)
model_direct1
summary(model_direct1)

# Step 2: Path A (X on M); SLC25A1 on PFKM
model_mediate1 <- lm(PFKM_counts ~ SLC25A1_counts, data=combo_df)
model_mediate1
summary(model_mediate1)

# Step 3: M on X, controlling for X; pathway ~ SLC25A1 + PFKM
# Finally, step3!
model_indirect1 <- lm(combo_df[,232] ~ SLC25A1_counts + PFKM_counts, data=combo_df)
model_indirect1
summary(model_indirect1)

# Is the mediation effect statistically significant?
citrate_results1 <- mediate(model_mediate1, model_indirect1, treat="SLC25A1_counts", mediator='PFKM_counts', boot=TRUE, sims=500)
summary(citrate_results1)

############################### Mediation analysis on individual genes ################################################
################################# Voltage dependent anion channel pathway ##############################################
# Find the activity scores for the VDAC pathway 
pathways <- colnames(act)
pathways
res <- str_detect(pathways, "tca")
which(res, arr.ind=TRUE) # 231
pathways[231]

# First, get the dataframe in the right shape
subdf <- horizontal_df %>% split(grepl('VDAC', x=.$gene)) # I can't find this gene for some reason...dropping it from mediation
head(subdf[[2]])

subdf1 <- horizontal_df %>% split(grepl('\\bSLC25A4\\b', x=.$gene)) # Check the subset dfs are the same nrow
head(subdf1[[2]])

subdf2 <- horizontal_df %>% split(grepl('\\bSLC25A5\\b', x=.$gene)) # Check the subset dfs are the same nrow
head(subdf2[[2]])

subdf3 <- horizontal_df %>% split(grepl('\\bSLC25A6\\b', x=.$gene)) # Check the subset dfs are the same nrow
head(subdf3[[2]])

subdf4 <- horizontal_df %>% split(grepl('\\bSLC25A31\\b', x=.$gene)) # Check the subset dfs are the same nrow
head(subdf4[[2]])

subdf5 <- horizontal_df %>% split(grepl('HK2', x=.$gene)) # Check the subset dfs are the same nrow
head(subdf5[[2]])

subdf6 <- horizontal_df %>% split(grepl('ND1', x=.$gene)) # Mitochondrial complex I --> I just chose one gene 
head(subdf6[[2]])

nrow(subdf[[2]])==nrow(subdf[[2]]) # TRUE

# Add this column to the subdf1 column
combo_df <- bind_rows(list(as.tibble(subdf1[[2]]), as.tibble(subdf2[[2]]), as.tibble(subdf3[[2]]), as.tibble(subdf4[[2]]), as.tibble(subdf5[[2]]), as.tibble(subdf6[[2]])))
combo_df[1:5,1:5]
tail(combo_df)
colnames(combo_df)
head(combo_df$gene)

# Rename the counts column in the first df so I can bind it to the second one
combo_df <- as.data.frame(combo_df)

colnames(combo_df)[296] <- "SLC25A4_counts" 
tail(combo_df)

colnames(combo_df)[299] <- "SLC25A5_counts" 
tail(combo_df)

colnames(combo_df)[300] <- "SLC25A6_counts" 
tail(combo_df)

colnames(combo_df)[301] <- "SLC25A31_counts" 
tail(combo_df)

colnames(combo_df)[302] <- "HK2_counts" 
tail(combo_df)

colnames(combo_df)[303] <- "ND1_counts" 
tail(combo_df)

colnames(combo_df)

# For some reason it put NA in the PFKM counts so I am going to force write it
combo_df$SLC25A4_counts <- subdf1[[2]]$SLC25A4_counts
combo_df$SLC25A5_counts <- subdf2[[2]]$SLC25A5_counts
#combo_df$SLC25A6_counts <- subdf3[[2]]$SLC25A6_counts # replacement has 592 rows, data has 2072
# dropping SLC25A6 as a mediatior for now...
combo_df$SLC25A31_counts <- subdf4[[2]]$SLC25A31_counts
combo_df$HK2_counts <- subdf5[[2]]$HK2_counts
combo_df$ND1_counts <- subdf6[[2]]$ND1_counts
tail(combo_df)

# Step 1: Estimate total effect (x on Y) ; HK2 on the TCA cycle   
model_direct1 <- lm(combo_df[,232] ~ SLC25A4_counts, data=combo_df)
model_direct1
summary(model_direct1)

# Step 2: Path A (X on M); complex I on HK2 
model_mediate1 <- lm(SLC25A4_counts ~ ND1_counts, data=combo_df)
model_mediate1
summary(model_mediate1)

# Step 3: M on X, controlling for X; pathway ~ SLC25A4 + HK2 
# Finally, step3!
model_indirect1 <- lm(combo_df[,232] ~ ND1_counts + SLC25A4_counts, data=combo_df)
model_indirect1
summary(model_indirect1)

# Is the mediation effect statistically significant?
citrate_results <- mediate(model_mediate1, model_indirect1, treat="ND1_counts", mediator='SLC25A4_counts', boot=TRUE, sims=500)
summary(citrate_results)


