# Script to do mediation analysis on the activity scores 

# Library
library(data.table)
library(janitor)
library(tidyverse)
#library(dplyr)
library(stringr)
library(mediation)
library(purrr)

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

# Mediation analysis
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
## Voltage dependent anion channel pathway
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
