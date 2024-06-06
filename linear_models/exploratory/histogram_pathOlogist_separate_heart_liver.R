# Script to do mediation analysis on the activity scores 

# Library
library(data.table)
library(janitor)
library(tidyverse)
library(dplyr)
library(stringr)
library(mediation)
library(purrr)
library(broom)
library(patchwork)
library(ggtext)

# Read in the pathways and activity scores at 50% activity
#activity <- fread("~/scratch/mjpete11/mitochondrial_paralogs/linear_models/data/data/filtered_activity_scores_50percent.csv",sep=",")
#activity[1:5,1:5]

# Read in the activity scores calculated for heart and liver separately
liver_activity <- read.table("/scratch/mjpete11/mitochondrial_paralogs/linear_models/data/data/Mollie_liver_split_pathway_activity.txt",sep="\t") 
heart_activity <- read.table("/scratch/mjpete11/mitochondrial_paralogs/linear_models/data/data/Mollie_heart_split_Pathway_Activity.txt",sep="\t") 

liver_activity[1:5,1:5]
heart_activity[1:5,1:5]

# Drop the first row and make the second row the colnames in each activity df
colnames(liver_activity)
colnames(liver_activity) <- NULL

colnames(heart_activity)
colnames(heart_activity) <- NULL

liver_activity <- liver_activity %>% row_to_names(row_number=1)
heart_activity <- heart_activity %>% row_to_names(row_number=1)

# Change the .s to   <-  in the sample names in the activity dataframe
#colnames(activity) <- gsub(x=colnames(activity), pattern="\\.", replacement="-")
#activity[1:5,1:5]

colnames(liver_activity) <- gsub(x=colnames(liver_activity), pattern="\\.", replacement="-")
liver_activity[1:5,1:5]

colnames(heart_activity) <- gsub(x=colnames(heart_activity), pattern="\\.", replacement="-")
heart_activity[1:5,1:5]

# rename the 'Pathway Name" column to snakecase
colnames(liver_activity)[1] <- c('pathway_name')
colnames(heart_activity)[1] <- c('pathway_name')

# Add a column with indicating the organ
liver_activity$organ <- "liver"
liver_activity <- liver_activity %>% dplyr::select(organ, everything())
dim(liver_activity) # 1324 150

heart_activity$organ <- "heart"
heart_activity <- heart_activity %>% dplyr::select(organ, everything())
dim(heart_activity) # 1324 150

# Rename the pathway column so it doesn't have a space
heart_activity$'Pathway Name' <- c('pathway_name')
heart_activity$pathway_name

liver_activity$'Pathway Name' <- c('pathway_name')
liver_activity$pathway_name

# Is there a perfect overlap in the pathway activities ; yes
length(intersect(heart_activity$pathway_name, liver_activity$pathway_name)) # 1324

# Combine the activity score dataframes
activity <- bind_rows(liver_activity, heart_activity) # this will vertical bind the dfs and fill in non-overlapping rows with NA
activity[1:5,1:5]
activity[1329:1324,1:5]
activity[1329:1324,151:156]
dim(activity) # 2648 299 

# Read in the gene counts with the RIN and activity scores
#organs <- read.csv("~/scratch/mjpete11/mitochondrial_paralogs/linear_models/data/data/organs_biomarkers_combat_seq3.csv", sep=",")%>%unique()
organs <- read.csv("/scratch/mjpete11/mitochondrial_paralogs/linear_models/data/data/qnorm_voom_SLC25_biomarkers_manifest.csv", sep=",")
organs[1:5,1:5]
dim(organs) # 34632 10

# List of SLC25 paralogs
SLC <- c("SLC25A1", "SLC25A2", "SLC25A3", "SLC25A4", "SLC25A5", "SLC25A6","UCP1", "UCP2", "UCP3", "SLC25A10", "SLC25A11", "SLC25A12","SLC25A13", "SLC25A14", "SLC25A15", "SLC25A16", "SLC25A17","SLC25A18", "SLC25A19", "SLC25A20", "SLC25A21", "SLC25A22","SLC25A23", "SLC25A24", "SLC25A25", "SLC25A26", "SLC25A27","SLC25A28", "SLC25A29", "SLC25A30", "SLC25A31", "SLC25A32","SLC25A33", "SLC25A34", "SLC25A35", "SLC25A36", "SLC25A37","SLC25A38", "SLC25A39", "SLC25A40", "SLC25A41", "SLC25A42","SLC25A43", "SLC25A44", "SLC25A45", "SLC25A46", "SLC25A47","SLC25A48", "MTCH1", "MTCH2", "SLC25A51", "SLC25A52", "SLC25A53")

# List of metabolic biomarkers 
biomarkers <- c("SLC16A1", "SLC16A2", "HK1", "HK2", "HK3", "GPI", "PFKM", "PFKP", "ALDOA", "ALDOB", "TPI1", "GAPDH", "PGK1","PGAM1", "PGAM5", "PGAM4", "ENO1","PKLR", "LDHA", "LDHB", "LDHC", "LDHD", "CPT1A","CPT1B", "CPT2", "ACADVL", "ACADL", "ACADM", "ACADS","ACADSB", "ACAD11", "ACAD8", "ACAD10", "ACAD9", "ECHS1","ECH1", "ECI1", "ECI2", "ECHS1", "HADHA", "PDHA1", "PDHA2","PDHB", "PDP1", "PDHX", "CS", "ACO1", "ACO2", "IDH3G","IDH1", "IDH2", "IDH3B", "OGDH", "SUCLG2","SUCLG1", "SUCLA2", "SDHB", "SDHD", "SDHA", "SDHC", "SDHAF2","FH", "MDH1", "MDH2", "MTND1")

# Combine list
genes <- c(SLC, biomarkers)
genes
length(genes) # 118

# Subset just the TCA cycle pathway activity scores since it is the dependent variable of interest
TCA_DF <- activity %>% filter(str_detect(activity$pathway_name,"tca"))
head(TCA_DF)
dim(TCA_DF) # 2 299
TCA_DF[1:5,1:5]
TCA_DF <- TCA_DF[c(1:2),]
TCA_DF[,1:5]

# Reshape df into tall format
TCA_DF_long <- melt(setDT(TCA_DF), id.vars=c("organ", "pathway_name"), variable.name="sample_ID")
head(TCA_DF_long)
tail(TCA_DF_long)

# Drop rows with NA (these are rows where the sample ID will not have an activity score because it does
# not match the organ column so melt() fills these values with NA)
TCA_DF_long <- TCA_DF_long %>% drop_na()

# Change the name of the value name to activity
colnames(TCA_DF_long)[4] <- c("activity") 

# Convert the activity column from character to numeric
TCA_DF_long$activity <- as.numeric(TCA_DF_long$activity)
class(TCA_DF_long$activity)
any(is.na(TCA_DF_long))==TRUE # TRUE
which(is.na(TCA_DF_long))
head(TCA_DF_long)
tail(TCA_DF_long)
TCA_DF_long <- TCA_DF_long %>% drop <- na()
any(is.na(TCA_DF_long))==TRUE # FALSE


# Make a dataframe to make a heatmap
HeatMAP <- organs %>% left_join(TCA_DF_long %>% dplyr::select(sample_ID, activity))

# Categorize if SLC25 or not
HeatMAP <- HeatMAP %>%
  mutate(GeneSet = if_else(str_detect(gene,"SLC"),"SLC25","Biomarker"))

# Remove the "-" in MT-DNA
HeatMAP <- HeatMAP %>% mutate(across(gene, str_replace, 'MT-ND1', 'MTND1'))
nrow(HeatMAP)==nrow(organs) # TRUE
tail(HeatMAP)

# Histograms of all the genes (SLC25 + biomarkers) after running PathOlogist on
# the heart and liver samples separately
# Histogram plot function
rm(plots)
rm(histogram_plots)
# Layered histogram plot
histogram_plots <- function(GENE){
  dat <- organs %>% filter(gene == GENE)
  # Violin plot
  #p <- ggplot(dat, aes(x = combat_seq_counts, fill=organ)) +
  p <- ggplot(dat, aes(x = qnorm_voom_logCPM, fill=organ)) +
    geom_histogram(bins=100) +
    scale_color_manual(values = c("red", "blue")) +
    #	geom_bin2d(bins = 100) +
    #ylim(c(0,30)) +
    ylab("frequency") +
    #	xlab("filtered and combat_seq adjusted counts") +
    xlab("filtered, quantile and voom adjusted counts") +
    theme(plot.title = element_textbox_simple(size=10)) + # automatically wrap long titles
    ggtitle(paste0("Histogram plot of ", GENE, " expression in heart and liver")) 
    ggsave(paste0("/scratch/mjpete11/mitochondrial_paralogs/linear_models/data/data/pathologist_separate_histogram_plots/", GENE, ".png"), p, device="png")
}
dev.new()
Map(histogram_plots, GENE=genes)
dev.off()

histogram_plots(GENE=genes[118])
