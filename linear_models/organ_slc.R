# Purpose: Linear models to asses if SLC paralogs are dependent on energetic state?

# Libraries
library(limma)
library(edgeR)
library(data.table)
library(tidyverse)
library(reshape)

# Read in GTEx manifest
manifest <- read.csv("sample.tsv", header=TRUE, sep = "\t")

# Make dataframe with sample id, tissue type
# All of the rin number and ischemic time values were missing...
df1 <- data.frame(manifest$"dbgap_sample_id", manifest$"tissue_type")

# Remove all rows with NA in either columns
df2 <- df1[complete.cases(df1), ]

# Read in GTEx counts
counts <- fread("GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct", sep="\t")

# List of SLC25 paralogs
SLC <- c("SLC25A1", "SLC25A2", "SLC25A3", "SLC25A4", "SLC25A5", "SLC25A6", 
		 "SLC25A7", "SLC25A8", "SLC25A9", "SLC25A10", "SLC25A11", "SLC25A12", 
		 "SLC25A13", "SLC25A14", "SLC25A15", "SLC25A16", "SLC25A17", 
		 "SLC25A18", "SLC25A19", "SLC25A20", "SLC25A21", "SLC25A22", 
		 "SLC25A23", "SLC25A24", "SLC25A25", "SLC25A26", "SLC25A27", 
		 "SLC25A28", "SLC25A29", "SLC25A30", "SLC25A31", "SLC25A32", 
		 "SLC25A33", "SLC25A34", "SLC25A35", "SLC25A36", "SLC25A37", 
		 "SLC25A38", "SLC25A39", "SLC25A40", "SLC25A41", "SLC25A42", 
		 "SLC25A43", "SLC25A44", "SLC25A45", "SLC25A46", "SLC25A47",
		 "SLC25A48", "SLC25A49", "SLC25A50", "SLC25A51", "SLC25A52", "SLC25A53")

# Subset the SLC25 genes
sub_df <- counts[counts$"Description" %in% SLC, ]

# Only 49 genes are present; which are missing?
setdiff(SLC, sub_df$'Description') 
# "SLC25A7"  "SLC25A8"  "SLC25A9"  "SLC25A49" "SLC25A50"
# I searched for the ENSG ID but I couldn't find them

# Read in sample attributes files
file2 <- read.csv("GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt", sep="\t")

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
heart4 <- melt(data = heart3,
			   id.vars = c("gene", "organ"),
			   measure.vars = colnames(heart3)[3:ncol(heart3)-1],
			   variable.name = "samples",
			   value.name = "counts")

liver4 <- melt(data = liver3,
			   id.vars = c("gene", "organ"),
			   measure.vars = colnames(liver3)[3:ncol(liver3)-1],
			   variable.name = "samples",
			   value.name = "counts")

# Combine into one df
organs <- rbind(heart4, liver4)

# Factor response variable
organs$organ <- as.factor(organs$organ)

# Linear model: All SLC ~ organ
model <- lm(value ~ 0 + organ, organs)
res <- summary(model)

# Write linear model results to file
#sink(paste("/scratch/mjpete11/linear_models/results/","SLC_organ.csv",sep=""))
#print(res)
#sink()

# Function to make linear model of each SLC ~ organ individually
# List of subset dfs
SLC_dfs <- c("A1", "A2", "A3", "A4", "A5", "A6", "A7", "A10",
			 "A11", "A12", "A13", "A14", "A15", "A16", "A17", "A18", "A19",
			 "A20", "A21", "A22", "A23", "A24", "A25", "A26", "A27", "A28",
			 "A29", "A30", "A31", "A32", "A33", "A34", "A35", "A36", "A37",
			 "A38", "A39", "A40", "A41", "A42", "A43", "A44", "A45", "A46",
			 "A47", "A48", "A51", "A52", "A53")

SLCs <- c("SLC25A1", "SLC25A2", "SLC25A3", "SLC25A4", "SLC25A5", "SLC25A6", 
		 "SLC25A7", "SLC25A10", "SLC25A11", "SLC25A12", 
		 "SLC25A13", "SLC25A14", "SLC25A15", "SLC25A16", "SLC25A17", 
		 "SLC25A18", "SLC25A19", "SLC25A20", "SLC25A21", "SLC25A22", 
		 "SLC25A23", "SLC25A24", "SLC25A25", "SLC25A26", "SLC25A27", 
		 "SLC25A28", "SLC25A29", "SLC25A30", "SLC25A31", "SLC25A32", 
		 "SLC25A33", "SLC25A34", "SLC25A35", "SLC25A36", "SLC25A37", 
		 "SLC25A38", "SLC25A39", "SLC25A40", "SLC25A41", "SLC25A42", 
		 "SLC25A43", "SLC25A44", "SLC25A45", "SLC25A46", "SLC25A47",
		 "SLC25A48", "SLC25A51", "SLC25A52", "SLC25A53")

models <- sprintf("model%d",1:49)
results <- sprintf("res%d",1:49)

# Function
#linear_model <- function(var1, var2, var3, var4, var5){
#	var1 <- organs %>% filter(gene == var2)
#	var3 <- lapply(var1, function(x) as.factor(x[["organ"]])) 
#	var4 <- lm(value ~ 0 + organ, var3)
#	var5 <- summary(var4)
#	sink(paste("/scratch/mjpete11/linear_models/results/", paste0(var2,"_organ.csv"),sep=""))
#	print(var5)
#	sink()
#}
#obj <- Map(linear_model,var1=SLC_dfs, var2=SLCs, var3=SLC_dfs) 
#obj <- Map(linear_model,var1=SLC_dfs, var2=SLCs, var3=models, var4=results) 

# Subset one gene from the organ df
A1 <- organs %>% filter(gene == "SLC25A1")
model1 <- lm(value ~ 0 + organ, A1)
res1 <- summary(model1)
sink(paste("/scratch/mjpete11/linear_models/results/","SLC1_organ.csv",sep=""))
print(res1)
sink()

# Subset one gene from the organ df
A2 <- organs %>% filter(gene == "SLC25A2")
model2 <- lm(value ~ 0 + organ, A2)
res2 <- summary(model2)
sink(paste("/scratch/mjpete11/linear_models/results/","SLC2_organ.csv",sep=""))
print(res2)
sink()

# Subset one gene from the organ df
A3 <- organs %>% filter(gene == "SLC25A3")
model3 <- lm(value ~ 0 + organ, A3)
res3 <- summary(model3)
sink(paste("/scratch/mjpete11/linear_models/results/","SLC3_organ.csv",sep=""))
print(res3)
sink()

# Subset one gene from the organ df
A4 <- organs %>% filter(gene == "SLC25A4")
model4 <- lm(value ~ 0 + organ, A4)
res4 <- summary(model4)
sink(paste("/scratch/mjpete11/linear_models/results/","SLC4_organ.csv",sep=""))
print(res4)
sink()

# Subset one gene from the organ df
A5 <- organs %>% filter(gene == "SLC25A5")
model5 <- lm(value ~ 0 + organ, A5)
res5 <- summary(model5)
sink(paste("/scratch/mjpete11/linear_models/results/","SLC5_organ.csv",sep=""))
print(res5)
sink()

# Subset one gene from the organ df
A6 <- organs %>% filter(gene == "SLC25A6")
model6 <- lm(value ~ 0 + organ, A6)
res6 <- summary(model6)
sink(paste("/scratch/mjpete11/linear_models/results/","SLC6_organ.csv",sep=""))
print(res6)
sink()

# Subset one gene from the organ df
A7 <- organs %>% filter(gene == "SLC25A10")
model7 <- lm(value ~ 0 + organ, A7)
res7 <- summary(model7)
sink(paste("/scratch/mjpete11/linear_models/results/","SLC10_organ.csv",sep=""))
print(res7)
sink()

# Subset one gene from the organ df
A8 <- organs %>% filter(gene == "SLC25A11")
model8 <- lm(value ~ 0 + organ, A8)
res8 <- summary(model8)
sink(paste("/scratch/mjpete11/linear_models/results/","SLC11_organ.csv",sep=""))
print(res8)
sink()

# Subset one gene from the organ df
A9 <- organs %>% filter(gene == "SLC25A12")
model9 <- lm(value ~ 0 + organ, A9)
res9 <- summary(model9)
sink(paste("/scratch/mjpete11/linear_models/results/","SLC12_organ.csv",sep=""))
print(res9)
sink()

# Subset one gene from the organ df
A10 <- organs %>% filter(gene == "SLC25A13")
model10 <- lm(value ~ 0 + organ, A10)
res10 <- summary(model10)
sink(paste("/scratch/mjpete11/linear_models/results/","SLC13_organ.csv",sep=""))
print(res10)
sink()

# Subset one gene from the organ df
A11 <- organs %>% filter(gene == "SLC25A14")
model11 <- lm(value ~ 0 + organ, A11)
res11 <- summary(model11)
sink(paste("/scratch/mjpete11/linear_models/results/","SLC14_organ.csv",sep=""))
print(res11)
sink()

# Subset one gene from the organ df
A12 <- organs %>% filter(gene == "SLC25A15")
model12 <- lm(value ~ 0 + organ, A12)
res12 <- summary(model12)
sink(paste("/scratch/mjpete11/linear_models/results/","SLC15_organ.csv",sep=""))
print(res12)
sink()

# Subset one gene from the organ df
A13 <- organs %>% filter(gene == "SLC25A16")
model13 <- lm(value ~ 0 + organ, A13)
res13 <- summary(model13)
sink(paste("/scratch/mjpete11/linear_models/results/","SLC16_organ.csv",sep=""))
print(res13)
sink()

# Subset one gene from the organ df
A14 <- organs %>% filter(gene == "SLC25A17")
model14 <- lm(value ~ 0 + organ, A14)
res14 <- summary(model14)
sink(paste("/scratch/mjpete11/linear_models/results/","SLC17_organ.csv",sep=""))
print(res14)
sink()

# Subset one gene from the organ df
A15 <- organs %>% filter(gene == "SLC25A18")
model15 <- lm(value ~ 0 + organ, A15)
res15 <- summary(model15)
sink(paste("/scratch/mjpete11/linear_models/results/","SLC18_organ.csv",sep=""))
print(res15)
sink()


# Subset one gene from the organ df
A16 <- organs %>% filter(gene == "SLC25A19")
model16 <- lm(value ~ 0 + organ, A16)
res16 <- summary(model16)
sink(paste("/scratch/mjpete11/linear_models/results/","SLC19_organ.csv",sep=""))
print(res16)
sink()

# Subset one gene from the organ df
A17 <- organs %>% filter(gene == "SLC25A20")
model17 <- lm(value ~ 0 + organ, A17)
res17 <- summary(model17)
sink(paste("/scratch/mjpete11/linear_models/results/","SLC20_organ.csv",sep=""))
print(res17)
sink()

# Subset one gene from the organ df
A18 <- organs %>% filter(gene == "SLC25A21")
model18 <- lm(value ~ 0 + organ, A18)
res18 <- summary(model18)
sink(paste("/scratch/mjpete11/linear_models/results/","SLC21_organ.csv",sep=""))
print(res18)
sink()

# Subset one gene from the organ df
A19 <- organs %>% filter(gene == "SLC25A22")
model19 <- lm(value ~ 0 + organ, A19)
res19 <- summary(model19)
sink(paste("/scratch/mjpete11/linear_models/results/","SLC22_organ.csv",sep=""))
print(res19)
sink()

# Subset one gene from the organ df
A20 <- organs %>% filter(gene == "SLC25A23")
model20 <- lm(value ~ 0 + organ, A20)
res20 <- summary(model20)
sink(paste("/scratch/mjpete11/linear_models/results/","SLC23_organ.csv",sep=""))
print(res20)
sink()

# Subset one gene from the organ df
A21 <- organs %>% filter(gene == "SLC25A24")
model21 <- lm(value ~ 0 + organ, A21)
res21 <- summary(model21)
sink(paste("/scratch/mjpete11/linear_models/results/","SLC24_organ.csv",sep=""))
print(res21)
sink()

# Subset one gene from the organ df
A22 <- organs %>% filter(gene == "SLC25A25")
model22 <- lm(value ~ 0 + organ, A22)
res22 <- summary(model22)
sink(paste("/scratch/mjpete11/linear_models/results/","SLC25_organ.csv",sep=""))
print(res22)
sink()

# Subset one gene from the organ df
A23 <- organs %>% filter(gene == "SLC25A26")
model23 <- lm(value ~ 0 + organ, A23)
res23 <- summary(model23)
sink(paste("/scratch/mjpete11/linear_models/results/","SLC26_organ.csv",sep=""))
print(res23)
sink()

# Subset one gene from the organ df
A24 <- organs %>% filter(gene == "SLC25A27")
model24 <- lm(value ~ 0 + organ, A24)
res24 <- summary(model24)
sink(paste("/scratch/mjpete11/linear_models/results/","SLC27_organ.csv",sep=""))
print(res24)
sink()

# Subset one gene from the organ df
A25 <- organs %>% filter(gene == "SLC25A28")
model25 <- lm(value ~ 0 + organ, A25)
res25 <- summary(model25)
sink(paste("/scratch/mjpete11/linear_models/results/","SLC28_organ.csv",sep=""))
print(res25)
sink()

# Subset one gene from the organ df
A26 <- organs %>% filter(gene == "SLC25A29")
model26 <- lm(value ~ 0 + organ, A26)
res26 <- summary(model26)
sink(paste("/scratch/mjpete11/linear_models/results/","SLC29_organ.csv",sep=""))
print(res26)
sink()

# Subset one gene from the organ df
A27 <- organs %>% filter(gene == "SLC25A30")
model27 <- lm(value ~ 0 + organ, A27)
res27 <- summary(model27)
sink(paste("/scratch/mjpete11/linear_models/results/","SLC30_organ.csv",sep=""))
print(res27)
sink()

# Subset one gene from the organ df
A28 <- organs %>% filter(gene == "SLC25A31")
model28 <- lm(value ~ 0 + organ, A28)
res28 <- summary(model28)
sink(paste("/scratch/mjpete11/linear_models/results/","SLC31_organ.csv",sep=""))
print(res28)
sink()

# Subset one gene from the organ df
A29 <- organs %>% filter(gene == "SLC25A32")
model29 <- lm(value ~ 0 + organ, A29)
res29 <- summary(model29)
sink(paste("/scratch/mjpete11/linear_models/results/","SLC32_organ.csv",sep=""))
print(res29)
sink()

# Subset one gene from the organ df
A30 <- organs %>% filter(gene == "SLC25A33")
model30 <- lm(value ~ 0 + organ, A30)
res30 <- summary(model30)
sink(paste("/scratch/mjpete11/linear_models/results/","SLC33_organ.csv",sep=""))
print(res30)
sink()

# Subset one gene from the organ df
A31 <- organs %>% filter(gene == "SLC25A34")
model31 <- lm(value ~ 0 + organ, A31)
res31 <- summary(model31)
sink(paste("/scratch/mjpete11/linear_models/results/","SLC34_organ.csv",sep=""))
print(res31)
sink()

# Subset one gene from the organ df
A32 <- organs %>% filter(gene == "SLC25A35")
model32 <- lm(value ~ 0 + organ, A32)
res32 <- summary(model32)
sink(paste("/scratch/mjpete11/linear_models/results/","SLC35_organ.csv",sep=""))
print(res32)
sink()

# Subset one gene from the organ df
A33 <- organs %>% filter(gene == "SLC25A36")
model33 <- lm(value ~ 0 + organ, A33)
res33 <- summary(model33)
sink(paste("/scratch/mjpete11/linear_models/results/","SLC36_organ.csv",sep=""))
print(res33)
sink()

# Subset one gene from the organ df
A34 <- organs %>% filter(gene == "SLC25A37")
model34 <- lm(value ~ 0 + organ, A34)
res34 <- summary(model34)
sink(paste("/scratch/mjpete11/linear_models/results/","SLC37_organ.csv",sep=""))
print(res34)
sink()

# Subset one gene from the organ df
A35 <- organs %>% filter(gene == "SLC25A38")
model35 <- lm(value ~ 0 + organ, A35)
res35 <- summary(model35)
sink(paste("/scratch/mjpete11/linear_models/results/","SLC38_organ.csv",sep=""))
print(res35)
sink()

# Subset one gene from the organ df
A36 <- organs %>% filter(gene == "SLC25A39")
model36 <- lm(value ~ 0 + organ, A36)
res36 <- summary(model36)
sink(paste("/scratch/mjpete11/linear_models/results/","SLC39_organ.csv",sep=""))
print(res36)
sink()

# Subset one gene from the organ df
A36 <- organs %>% filter(gene == "SLC25A39")
model36 <- lm(value ~ 0 + organ, A36)
res36 <- summary(model36)
sink(paste("/scratch/mjpete11/linear_models/results/","SLC39_organ.csv",sep=""))
print(res36)
sink()

# Subset one gene from the organ df
A37 <- organs %>% filter(gene == "SLC25A40")
model37 <- lm(value ~ 0 + organ, A37)
res37 <- summary(model37)
sink(paste("/scratch/mjpete11/linear_models/results/","SLC40_organ.csv",sep=""))
print(res37)
sink()

# Subset one gene from the organ df
A38 <- organs %>% filter(gene == "SLC25A41")
model38 <- lm(value ~ 0 + organ, A38)
res38 <- summary(model38)
sink(paste("/scratch/mjpete11/linear_models/results/","SLC41_organ.csv",sep=""))
print(res38)
sink()

# Subset one gene from the organ df
A42 <- organs %>% filter(gene == "SLC25A42")
model42 <- lm(value ~ 0 + organ, A42)
res42 <- summary(model42)
sink(paste("/scratch/mjpete11/linear_models/results/","SLC42_organ.csv",sep=""))
print(res42)
sink()

# Subset one gene from the organ df
A43 <- organs %>% filter(gene == "SLC25A43")
model43 <- lm(value ~ 0 + organ, A43)
res43 <- summary(model43)
sink(paste("/scratch/mjpete11/linear_models/results/","SLC43_organ.csv",sep=""))
print(res43)
sink()

# Subset one gene from the organ df
A44 <- organs %>% filter(gene == "SLC25A44")
model44 <- lm(value ~ 0 + organ, A44)
res44 <- summary(model44)
sink(paste("/scratch/mjpete11/linear_models/results/","SLC44_organ.csv",sep=""))
print(res44)
sink()

# Subset one gene from the organ df
A45 <- organs %>% filter(gene == "SLC25A45")
model45 <- lm(value ~ 0 + organ, A45)
res45 <- summary(model45)
sink(paste("/scratch/mjpete11/linear_models/results/","SLC45_organ.csv",sep=""))
print(res45)
sink()

# Subset one gene from the organ df
A46 <- organs %>% filter(gene == "SLC25A46")
model46 <- lm(value ~ 0 + organ, A46)
res46 <- summary(model46)
sink(paste("/scratch/mjpete11/linear_models/results/","SLC46_organ.csv",sep=""))
print(res46)
sink()

# Subset one gene from the organ df
A47 <- organs %>% filter(gene == "SLC25A47")
model47 <- lm(value ~ 0 + organ, A47)
res47 <- summary(model47)
sink(paste("/scratch/mjpete11/linear_models/results/","SLC47_organ.csv",sep=""))
print(res47)
sink()

# Subset one gene from the organ df
A48 <- organs %>% filter(gene == "SLC25A48")
model48 <- lm(value ~ 0 + organ, A48)
res48 <- summary(model48)
sink(paste("/scratch/mjpete11/linear_models/results/","SLC48_organ.csv",sep=""))
print(res48)
sink()

# Subset one gene from the organ df
A49 <- organs %>% filter(gene == "SLC25A51")
model49 <- lm(value ~ 0 + organ, A49)
res51 <- summary(model49)
sink(paste("/scratch/mjpete11/linear_models/results/","SLC51_organ.csv",sep=""))
print(res51)
sink()

# Subset one gene from the organ df
A52 <- organs %>% filter(gene == "SLC25A52")
model52 <- lm(value ~ 0 + organ, A52)
res52 <- summary(model52)
sink(paste("/scratch/mjpete11/linear_models/results/","SLC52_organ.csv",sep=""))
print(res52)
sink()

# Subset one gene from the organ df
A53 <- organs %>% filter(gene == "SLC25A53")
model53 <- lm(value ~ 0 + organ, A53)
res53 <- summary(model53)
sink(paste("/scratch/mjpete11/linear_models/results/","SLC53_organ.csv",sep=""))
print(res53)
sink()
