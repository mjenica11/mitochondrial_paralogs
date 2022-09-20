# Purpose: Linear models to asses if SLC paralogs are dependent on energetic state?

# Libraries
library(limma)
library(edgeR)
library(data.table)
library(tidyverse)
library(reshape)
library(plyr)
library(moderndive)
library(infer)
library(performance)
library(see)
library(ggplot2)
library(ggpubr)

# Read in GTEx manifest
manifest <- read.csv("/scratch/mjpete11/linear_models/data/sample.tsv", header=TRUE, sep = "\t")

# Make dataframe with sample id, tissue type
# All of the rin number and ischemic time values were missing...
df1 <- data.frame(manifest$"dbgap_sample_id", manifest$"tissue_type")

# Remove all rows with NA in either columns
df2 <- df1[complete.cases(df1), ]

# Read in GTEx counts
counts <- fread("/scratch/mjpete11/linear_models/data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct", sep="\t")

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

# Change the name of the 'variable' column to 'SAMPID' to match columns
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
#organs <- rbind(heart4, liver4)
#organs <- left_join(heart4, liver4, by="gene")

# Add a column with the individual IDs
# Add sex and age to manifest
# Read in separate file
meta <- read.csv("/scratch/mjpete11/linear_models/data/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt", header=TRUE, sep = "\t")

# Add sample ID col to meta
library(stringr)
heart4$SUBJID <- str_sub(heart4$SAMPID, start=1L, end=10L)
liver4$SUBJID <- str_sub(liver4$SAMPID, start=1L, end=10L)

# Subset to only heart and liver samples that are from the same person
heart4$SMGEBTCH <- NULL
heart4$SMRIN <- NULL
heart4$SMTSISCH <- NULL
heart4$SMGEBTCHT <- NULL

liver4$SMGEBTCH <- NULL
liver4$SMRIN <- NULL
liver4$SMTSISCH <- NULL
liver4$SMGEBTCHT <- NULL

# Dataframe with only samples from people who donated both a heart and liver
organs <- merge(heart4, liver4, by="SUBJID") # 148 samples
tmp1 <- organs[,c("SUBJID", "gene.x", "organ.x", "SAMPID.x", "value.x")]
tmp2 <- organs[,c("SUBJID", "gene.y", "organ.y", "SAMPID.y", "value.y")]
colnames(tmp1) <- c("SUBJID", "gene", "organ", "SAMPID", "value")
colnames(tmp2) <- c("SUBJID", "gene", "organ", "SAMPID", "value")
organs2 <- rbind(tmp1, tmp2)

# Dataframe of all heart and liver samples
#organs <- rbind.fill(heart4, liver4)
#head(organs)

#Tab <- organs %>% select(SUBJID,organ)%>%unique()%>%pull(SUBJID)%>%table()
#both_tissue<-names(Tab)[which(as.numeric(Tab)==2)]
#test <- test %>% slice(which(SUBJID%in%both_tissue))
#head(test)
#tail(test)

# Reshape organs df 
#organs$SMGEBTCH <- NULL
#organs$SMRIN <- NULL
#organs$SMTSISCH <- NULL
#organs$SMGEBTCHT <- NULL

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

# Function
# Function to make box plots of each SLC ~ organ individually
# List of subset dfs
SLC <- dfs <- c("A1", "A2", "A3", "A4", "A5", "A6", "A7", "A10",
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
SLCs <- as.list(SLCs)
box_plots <- function(GENE){
			dat <- organs2 %>% filter(gene == GENE)
			dat <- na.omit(dat)
			p <- ggboxplot(dat, x = "organ", y = "value",
					   	   color = "organ", palette = c("#FFA500", "#FF0000"),
						   ylab = "counts", xlab = "organ")
		    p + stat_compare_means()
			ggsave(paste0(GENE, "_boxplot.pdf"))
}
#plots <- Map(box_plots, GENE=SLCs)


# Subset one gene from the organ df
A1 <- organs2 %>% filter(gene == "SLC25A1")
A1 <- na.omit(A1)
model1 <- lm(value ~ 0 + organ + SUBJID, data=A1)
res1 <- summary(model1)
sink(paste("/scratch/mjpete11/linear_models/results3/","SLC1_organ.csv",sep=""))
print(res1)
sink()

# Diagnostic plots 
m <- check_model(model1, show_dots=TRUE, se=FALSE)
check_normality(model1)
check_collinearity(model1)
pdf(paste("/scratch/mjpete11/linear_models/results3/","SLC1_check_model.pdf",sep=""))
m
dev.off()
					      
# Subset one gene from the organ df
A2 <- organs2 %>% filter(gene == "SLC25A2")
model2 <- lm(value ~ 0 + organ + SUBJID, A2)
res2 <- summary(model2)
sink(paste("/scratch/mjpete11/linear_models/results3/","SLC2_organ.csv",sep=""))
print(res2)
sink()

# Diagnostic plots 
m2 <- check_model(model2, show_dots=TRUE, se=FALSE)
check_normality(model2)
check_collinearity(model2)
pdf(paste("/scratch/mjpete11/linear_models/results3/","SLC2_check_model.pdf",sep=""))
m2
dev.off()

# Subset one gene from the organ df
A3 <- organs2 %>% filter(gene == "SLC25A3")
model3 <- lm(value ~ 0 + organ + SUBJID, A3)
res3 <- summary(model3)
sink(paste("/scratch/mjpete11/linear_models/results3/","SLC3_organ.csv",sep=""))
print(res3)
sink()

# Diagnostic plots 
m3 <- check_model(model3, show_dots=TRUE, se=FALSE)
check_normality(model3)
check_collinearity(model3)
pdf(paste("/scratch/mjpete11/linear_models/results3/","SLC3_check_model.pdf",sep=""))
m3
dev.off()

# Subset one gene from the organ df
A4 <- organs2 %>% filter(gene == "SLC25A4")
model4 <- lm(value ~ 0 + organ + SUBJID, A4)
res4 <- summary(model4)
sink(paste("/scratch/mjpete11/linear_models/results3/","SLC4_organ.csv",sep=""))
print(res4)
sink()

# Diagnostic plots 
m4 <- check_model(model4, show_dots=TRUE, se=FALSE)
check_normality(model4)
check_collinearity(model4)
pdf(paste("/scratch/mjpete11/linear_models/results3/","SLC4_check_model.pdf",sep=""))
m4
dev.off()

# Subset one gene from the organ df
A5 <- organs2 %>% filter(gene == "SLC25A5")
model5 <- lm(value ~ 0 + organ + SUBJID, A5)
res5 <- summary(model5)
sink(paste("/scratch/mjpete11/linear_models/results3/","SLC5_organ.csv",sep=""))
print(res5)
sink()

# Diagnostic plots 
m5 <- check_model(model5, show_dots=TRUE, se=FALSE)
check_normality(model5)
check_collinearity(model5)
pdf(paste("/scratch/mjpete11/linear_models/results3/","SLC5_check_model.pdf",sep=""))
m5
dev.off()

# Subset one gene from the organ df
A6 <- organs2 %>% filter(gene == "SLC25A6")
model6 <- lm(value ~ 0 + organ + SUBJID, A6)
res6 <- summary(model6)
sink(paste("/scratch/mjpete11/linear_models/results3/","SLC6_organ.csv",sep=""))
print(res6)
sink()

# Diagnostic plots 
m6 <- check_model(model6, show_dots=TRUE, se=FALSE)
check_normality(model6)
check_collinearity(model6)
pdf(paste("/scratch/mjpete11/linear_models/results3/","SLC6_check_model.pdf",sep=""))
m6
dev.off()

# Subset one gene from the organ df
A7 <- organs2 %>% filter(gene == "SLC25A10")
model7 <- lm(value ~ 0 + organ + SUBJID, A7)
res7 <- summary(model7)
sink(paste("/scratch/mjpete11/linear_models/results3/","SLC10_organ.csv",sep=""))
print(res7)
sink()

# Diagnostic plots 
m7 <- check_model(model7, show_dots=TRUE, se=FALSE)
check_normality(model7)
check_collinearity(model7)
pdf(paste("/scratch/mjpete11/linear_models/results3/","SLC10_check_model.pdf",sep=""))
m7
dev.off()

# Subset one gene from the organ df
A8 <- organs2 %>% filter(gene == "SLC25A11")
model8 <- lm(value ~ 0 + organ + SUBJID, A8)
res8 <- summary(model8)
sink(paste("/scratch/mjpete11/linear_models/results3/","SLC11_organ.csv",sep=""))
print(res8)
sink()

# Diagnostic plots 
m8 <- check_model(model8, show_dots=TRUE, se=FALSE)
check_normality(model8)
check_collinearity(model8)
pdf(paste("/scratch/mjpete11/linear_models/results3/","SLC11_check_model.pdf",sep=""))
m8
dev.off()

# Subset one gene from the organ df
A9 <- organs2 %>% filter(gene == "SLC25A12")
model9 <- lm(value ~ 0 + organ + SUBJID, A9)
res9 <- summary(model9)
sink(paste("/scratch/mjpete11/linear_models/results3/","SLC12_organ.csv",sep=""))
print(res9)
sink()

# Diagnostic plots 
m9 <- check_model(model9, show_dots=TRUE, se=FALSE)
check_normality(model9)
check_collinearity(model9)
pdf(paste("/scratch/mjpete11/linear_models/results3/","SLC12_check_model.pdf",sep=""))
m9
dev.off()

# Subset one gene from the organ df
A10 <- organs2 %>% filter(gene == "SLC25A13")
model10 <- lm(value ~ 0 + organ + SUBJID, A10)
res10 <- summary(model10)
sink(paste("/scratch/mjpete11/linear_models/results3/","SLC13_organ.csv",sep=""))
print(res10)
sink()

# Diagnostic plots 
m10 <- check_model(model10, show_dots=TRUE, se=FALSE)
check_normality(model10)
check_collinearity(model10)
pdf(paste("/scratch/mjpete11/linear_models/results3/","SLC13_check_model.pdf",sep=""))
m10
dev.off()

# Subset one gene from the organ df
A11 <- organs2 %>% filter(gene == "SLC25A14")
model11 <- lm(value ~ 0 + organ + SUBJID, A11)
res11 <- summary(model11)
sink(paste("/scratch/mjpete11/linear_models/results3/","SLC14_organ.csv",sep=""))
print(res11)
sink()

# Diagnostic plots 
m11 <- check_model(model11, show_dots=TRUE, se=FALSE)
check_normality(model11)
check_collinearity(model11)
pdf(paste("/scratch/mjpete11/linear_models/results3/","SLC14_check_model.pdf",sep=""))
m11
dev.off()

# Subset one gene from the organ df
A12 <- organs2 %>% filter(gene == "SLC25A15")
model12 <- lm(value ~ 0 + organ + SUBJID, A12)
res12 <- summary(model12)
sink(paste("/scratch/mjpete11/linear_models/results3/","SLC15_organ.csv",sep=""))
print(res12)
sink()

# Diagnostic plots 
m12 <- check_model(model12, show_dots=TRUE, se=FALSE)
check_normality(model12)
check_collinearity(model12)
pdf(paste("/scratch/mjpete11/linear_models/results3/","SLC15_check_model.pdf",sep=""))
m12
dev.off()

# Subset one gene from the organ df
A13 <- organs2 %>% filter(gene == "SLC25A16")
model13 <- lm(value ~ 0 + organ + SUBJID, A13)
res13 <- summary(model13)
sink(paste("/scratch/mjpete11/linear_models/results3/","SLC16_organ.csv",sep=""))
print(res13)
sink()

# Diagnostic plots 
m13 <- check_model(model13, show_dots=TRUE, se=FALSE)
check_normality(model13)
check_collinearity(model13)
pdf(paste("/scratch/mjpete11/linear_models/results3/","SLC16_check_model.pdf",sep=""))
m13
dev.off()

# Subset one gene from the organ df
A14 <- organs2 %>% filter(gene == "SLC25A17")
model14 <- lm(value ~ 0 + organ + SUBJID, A14)
res14 <- summary(model14)
sink(paste("/scratch/mjpete11/linear_models/results3/","SLC17_organ.csv",sep=""))
print(res14)
sink()

# Diagnostic plots 
m14 <- check_model(model14, show_dots=TRUE, se=FALSE)
check_normality(model14)
check_collinearity(model14)
pdf(paste("/scratch/mjpete11/linear_models/results3/","SLC17_check_model.pdf",sep=""))
m14
dev.off()

# Subset one gene from the organ df
A15 <- organs2 %>% filter(gene == "SLC25A18")
model15 <- lm(value ~ 0 + organ + SUBJID, A15)
res15 <- summary(model15)
sink(paste("/scratch/mjpete11/linear_models/results3/","SLC18_organ.csv",sep=""))
print(res15)
sink()

# Diagnostic plots 
m15 <- check_model(model15, show_dots=TRUE, se=FALSE)
check_normality(model15)
check_collinearity(model15)
pdf(paste("/scratch/mjpete11/linear_models/results3/","SLC18_check_model.pdf",sep=""))
m15
dev.off()

# Subset one gene from the organ df
A16 <- organs2 %>% filter(gene == "SLC25A19")
model16 <- lm(value ~ 0 + organ + SUBJID, A16)
res16 <- summary(model16)
sink(paste("/scratch/mjpete11/linear_models/results3/","SLC19_organ.csv",sep=""))
print(res16)
sink()

# Diagnostic plots 
m16 <- check_model(model16, show_dots=TRUE, se=FALSE)
check_normality(model16)
check_collinearity(model16)
pdf(paste("/scratch/mjpete11/linear_models/results3/","SLC19_check_model.pdf",sep=""))
m16
dev.off()

# Subset one gene from the organ df
A17 <- organs2 %>% filter(gene == "SLC25A20")
model17 <- lm(value ~ 0 + organ + SUBJID, A17)
res17 <- summary(model17)
sink(paste("/scratch/mjpete11/linear_models/results3/","SLC20_organ.csv",sep=""))
print(res17)
sink()

# Diagnostic plots 
m17 <- check_model(model17, show_dots=TRUE, se=FALSE)
check_normality(model17)
check_collinearity(model17)
pdf(paste("/scratch/mjpete11/linear_models/results3/","SLC20_check_model.pdf",sep=""))
m17
dev.off()

# Subset one gene from the organ df
A18 <- organs2 %>% filter(gene == "SLC25A21")
model18 <- lm(value ~ 0 + organ + SUBJID, A18)
res18 <- summary(model18)
sink(paste("/scratch/mjpete11/linear_models/results3/","SLC21_organ.csv",sep=""))
print(res18)
sink()

# Diagnostic plots 
m18 <- check_model(model18, show_dots=TRUE, se=FALSE)
check_normality(model18)
check_collinearity(model18)
pdf(paste("/scratch/mjpete11/linear_models/results3/","SLC21_check_model.pdf",sep=""))
m18
dev.off()

# Subset one gene from the organ df
A19 <- organs2 %>% filter(gene == "SLC25A22")
model19 <- lm(value ~ 0 + organ + SUBJID, A19)
res19 <- summary(model19)
sink(paste("/scratch/mjpete11/linear_models/results3/","SLC22_organ.csv",sep=""))
print(res19)
sink()

# Diagnostic plots 
m19 <- check_model(model19, show_dots=TRUE, se=FALSE)
check_normality(model19)
check_collinearity(model19)
pdf(paste("/scratch/mjpete11/linear_models/results3/","SLC22_check_model.pdf",sep=""))
m19
dev.off()

# Subset one gene from the organ df
A20 <- organs2 %>% filter(gene == "SLC25A23")
model20 <- lm(value ~ 0 + organ + SUBJID, A20)
res20 <- summary(model20)
sink(paste("/scratch/mjpete11/linear_models/results3/","SLC23_organ.csv",sep=""))
print(res20)
sink()

# Diagnostic plots 
m20 <- check_model(model20, show_dots=TRUE, se=FALSE)
check_normality(model20)
check_collinearity(model20)
pdf(paste("/scratch/mjpete11/linear_models/results3/","SLC23_check_model.pdf",sep=""))
m20
dev.off()

# Subset one gene from the organ df
A21 <- organs2 %>% filter(gene == "SLC25A24")
model21 <- lm(value ~ 0 + organ + SUBJID, A21)
res21 <- summary(model21)
sink(paste("/scratch/mjpete11/linear_models/results3/","SLC24_organ.csv",sep=""))
print(res21)
sink()

# Diagnostic plots 
m21 <- check_model(model21, show_dots=TRUE, se=FALSE)
check_normality(model21)
check_collinearity(model21)
pdf(paste("/scratch/mjpete11/linear_models/results3/","SLC24_check_model.pdf",sep=""))
m21
dev.off()

# Subset one gene from the organ df
A22 <- organs2 %>% filter(gene == "SLC25A25")
model22 <- lm(value ~ 0 + organ + SUBJID, A22)
res22 <- summary(model22)
sink(paste("/scratch/mjpete11/linear_models/results3/","SLC25_organ.csv",sep=""))
print(res22)
sink()

# Diagnostic plots 
m22 <- check_model(model22, show_dots=TRUE, se=FALSE)
check_normality(model22)
check_collinearity(model22)
pdf(paste("/scratch/mjpete11/linear_models/results3/","SLC25_check_model.pdf",sep=""))
m22
dev.off()

# Subset one gene from the organ df
A23 <- organs2 %>% filter(gene == "SLC25A26")
model23 <- lm(value ~ 0 + organ + SUBJID, A23)
res23 <- summary(model23)
sink(paste("/scratch/mjpete11/linear_models/results3/","SLC26_organ.csv",sep=""))
print(res23)
sink()

# Diagnostic plots 
m23 <- check_model(model23, show_dots=TRUE, se=FALSE)
check_normality(model23)
check_collinearity(model23)
pdf(paste("/scratch/mjpete11/linear_models/results3/","SLC26_check_model.pdf",sep=""))
m23
dev.off()

# Subset one gene from the organ df
A24 <- organs2 %>% filter(gene == "SLC25A27")
model24 <- lm(value ~ 0 + organ + SUBJID, A24)
res24 <- summary(model24)
sink(paste("/scratch/mjpete11/linear_models/results3/","SLC27_organ.csv",sep=""))
print(res24)
sink()

# Diagnostic plots 
m24 <- check_model(model24, show_dots=TRUE, se=FALSE)
check_normality(model24)
check_collinearity(model24)
pdf(paste("/scratch/mjpete11/linear_models/results3/","SLC27_check_model.pdf",sep=""))
m24
dev.off()

# Subset one gene from the organ df
A25 <- organs2 %>% filter(gene == "SLC25A28")
model25 <- lm(value ~ 0 + organ + SUBJID, A25)
res25 <- summary(model25)
sink(paste("/scratch/mjpete11/linear_models/results3/","SLC28_organ.csv",sep=""))
print(res25)
sink()

# Diagnostic plots 
m25 <- check_model(model25, show_dots=TRUE, se=FALSE)
check_normality(model25)
check_collinearity(model25)
pdf(paste("/scratch/mjpete11/linear_models/results3/","SLC28_check_model.pdf",sep=""))
m25
dev.off()

# Subset one gene from the organ df
A26 <- organs2 %>% filter(gene == "SLC25A29")
model26 <- lm(value ~ 0 + organ + SUBJID, A26)
res26 <- summary(model26)
sink(paste("/scratch/mjpete11/linear_models/results3/","SLC29_organ.csv",sep=""))
print(res26)
sink()

# Diagnostic plots 
m26 <- check_model(model26, show_dots=TRUE, se=FALSE)
check_normality(model26)
check_collinearity(model26)
pdf(paste("/scratch/mjpete11/linear_models/results3/","SLC29_check_model.pdf",sep=""))
m26
dev.off()

# Subset one gene from the organ df
A27 <- organs2 %>% filter(gene == "SLC25A30")
model27 <- lm(value ~ 0 + organ + SUBJID, A27)
res27 <- summary(model27)
sink(paste("/scratch/mjpete11/linear_models/results3/","SLC30_organ.csv",sep=""))
print(res27)
sink()

# Diagnostic plots 
m27 <- check_model(model27, show_dots=TRUE, se=FALSE)
check_normality(model27)
check_collinearity(model27)
pdf(paste("/scratch/mjpete11/linear_models/results3/","SLC30_check_model.pdf",sep=""))
m27
dev.off()

# Subset one gene from the organ df
A28 <- organs2 %>% filter(gene == "SLC25A31")
model28 <- lm(value ~ 0 + organ + SUBJID, A28)
res28 <- summary(model28)
sink(paste("/scratch/mjpete11/linear_models/results3/","SLC31_organ.csv",sep=""))
print(res28)
sink()

# Diagnostic plots 
m28 <- check_model(model28, show_dots=TRUE, se=FALSE)
check_normality(model28)
check_collinearity(model28)
pdf(paste("/scratch/mjpete11/linear_models/results3/","SLC31_check_model.pdf",sep=""))
m28
dev.off()

# Subset one gene from the organ df
A29 <- organs2 %>% filter(gene == "SLC25A32")
model29 <- lm(value ~ 0 + organ + SUBJID, A29)
res29 <- summary(model29)
sink(paste("/scratch/mjpete11/linear_models/results3/","SLC32_organ.csv",sep=""))
print(res29)
sink()

# Diagnostic plots 
m29 <- check_model(model29, show_dots=TRUE, se=FALSE)
check_normality(model29)
check_collinearity(model29)
pdf(paste("/scratch/mjpete11/linear_models/results3/","SLC32_check_model.pdf",sep=""))
m29
dev.off()

# Subset one gene from the organ df
A30 <- organs2 %>% filter(gene == "SLC25A33")
model30 <- lm(value ~ 0 + organ + SUBJID, A30)
res30 <- summary(model30)
sink(paste("/scratch/mjpete11/linear_models/results3/","SLC33_organ.csv",sep=""))
print(res30)
sink()

# Diagnostic plots 
m30 <- check_model(model30, show_dots=TRUE, se=FALSE)
check_normality(model30)
check_collinearity(model30)
pdf(paste("/scratch/mjpete11/linear_models/results3/","SLC33_check_model.pdf",sep=""))
m30
dev.off()

# Subset one gene from the organ df
A31 <- organs2 %>% filter(gene == "SLC25A34")
model31 <- lm(value ~ 0 + organ + SUBJID, A31)
res31 <- summary(model31)
sink(paste("/scratch/mjpete11/linear_models/results3/","SLC34_organ.csv",sep=""))
print(res31)
sink()

# Diagnostic plots 
m31 <- check_model(model31, show_dots=TRUE, se=FALSE)
check_normality(model31)
check_collinearity(model31)
pdf(paste("/scratch/mjpete11/linear_models/results3/","SLC34_check_model.pdf",sep=""))
m31
dev.off()

# Subset one gene from the organ df
A32 <- organs2 %>% filter(gene == "SLC25A35")
model32 <- lm(value ~ 0 + organ + SUBJID, A32)
res32 <- summary(model32)
sink(paste("/scratch/mjpete11/linear_models/results3/","SLC35_organ.csv",sep=""))
print(res32)
sink()

# Diagnostic plots 
m32 <- check_model(model32, show_dots=TRUE, se=FALSE)
check_normality(model32)
check_collinearity(model32)
pdf(paste("/scratch/mjpete11/linear_models/results3/","SLC35_check_model.pdf",sep=""))
m32
dev.off()

# Subset one gene from the organ df
A33 <- organs2 %>% filter(gene == "SLC25A36")
model33 <- lm(value ~ 0 + organ + SUBJID, A33)
res33 <- summary(model33)
sink(paste("/scratch/mjpete11/linear_models/results3/","SLC36_organ.csv",sep=""))
print(res33)
sink()

# Diagnostic plots 
m33 <- check_model(model33, show_dots=TRUE, se=FALSE)
check_normality(model33)
check_collinearity(model33)
pdf(paste("/scratch/mjpete11/linear_models/results3/","SLC36_check_model.pdf",sep=""))
m33
dev.off()

# Subset one gene from the organ df
A34 <- organs2 %>% filter(gene == "SLC25A37")
model34 <- lm(value ~ 0 + organ + SUBJID, A34)
res34 <- summary(model34)
sink(paste("/scratch/mjpete11/linear_models/results3/","SLC37_organ.csv",sep=""))
print(res34)
sink()

# Diagnostic plots 
m34 <- check_model(model19, show_dots=TRUE, se=FALSE)
check_normality(model34)
check_collinearity(model34)
pdf(paste("/scratch/mjpete11/linear_models/results3/","SLC37_check_model.pdf",sep=""))
m34
dev.off()

# Subset one gene from the organ df
A35 <- organs2 %>% filter(gene == "SLC25A38")
model35 <- lm(value ~ 0 + organ + SUBJID, A35)
res35 <- summary(model35)
sink(paste("/scratch/mjpete11/linear_models/results3/","SLC38_organ.csv",sep=""))
print(res35)
sink()

# Diagnostic plots 
m35 <- check_model(model35, show_dots=TRUE, se=FALSE)
check_normality(model35)
check_collinearity(model35)
pdf(paste("/scratch/mjpete11/linear_models/results3/","SLC38_check_model.pdf",sep=""))
m35
dev.off()

# Subset one gene from the organ df
A36 <- organs2 %>% filter(gene == "SLC25A39")
model36 <- lm(value ~ 0 + organ + SUBJID, A36)
res36 <- summary(model36)
sink(paste("/scratch/mjpete11/linear_models/results3/","SLC39_organ.csv",sep=""))
print(res36)
sink()

# Diagnostic plots 
m36 <- check_model(model36, show_dots=TRUE, se=FALSE)
check_normality(model36)
check_collinearity(model36)
pdf(paste("/scratch/mjpete11/linear_models/results3/","SLC39_check_model.pdf",sep=""))
m36
dev.off()

# Subset one gene from the organ df
A37 <- organs2 %>% filter(gene == "SLC25A40")
model37 <- lm(value ~ 0 + organ + SUBJID, A37)
res37 <- summary(model37)
sink(paste("/scratch/mjpete11/linear_models/results3/","SLC40_organ.csv",sep=""))
print(res37)
sink()

# Diagnostic plots 
m37 <- check_model(model37, show_dots=TRUE, se=FALSE)
check_normality(model37)
check_collinearity(model37)
pdf(paste("/scratch/mjpete11/linear_models/results3/","SLC40_check_model.pdf",sep=""))
m37
dev.off()

# Subset one gene from the organ df
A38 <- organs2 %>% filter(gene == "SLC25A41")
model38 <- lm(value ~ 0 + organ + SUBJID, A38)
res38 <- summary(model38)
sink(paste("/scratch/mjpete11/linear_models/results3/","SLC41_organ.csv",sep=""))
print(res38)
sink()

# Diagnostic plots 
m38 <- check_model(model38, show_dots=TRUE, se=FALSE)
check_normality(model38)
check_collinearity(model38)
pdf(paste("/scratch/mjpete11/linear_models/results3/","SLC41_check_model.pdf",sep=""))
m38
dev.off()

# Subset one gene from the organ df
A42 <- organs2 %>% filter(gene == "SLC25A42")
model42 <- lm(value ~ 0 + organ + SUBJID, A42)
res42 <- summary(model42)
sink(paste("/scratch/mjpete11/linear_models/results3/","SLC42_organ.csv",sep=""))
print(res42)
sink()

# Diagnostic plots 
m42 <- check_model(model42, show_dots=TRUE, se=FALSE)
check_normality(model42)
check_collinearity(model42)
pdf(paste("/scratch/mjpete11/linear_models/results3/","SLC42_check_model.pdf",sep=""))
m42
dev.off()

# Subset one gene from the organ df
A43 <- organs2 %>% filter(gene == "SLC25A43")
model43 <- lm(value ~ 0 + organ + SUBJID, A43)
res43 <- summary(model43)
sink(paste("/scratch/mjpete11/linear_models/results3/","SLC43_organ.csv",sep=""))
print(res43)
sink()

# Diagnostic plots 
m43 <- check_model(model43, show_dots=TRUE, se=FALSE)
check_normality(model43)
check_collinearity(model43)
pdf(paste("/scratch/mjpete11/linear_models/results3/","SLC43_check_model.pdf",sep=""))
m43
dev.off()

# Subset one gene from the organ df
A44 <- organs2 %>% filter(gene == "SLC25A44")
model44 <- lm(value ~ 0 + organ + SUBJID, A44)
res44 <- summary(model44)
sink(paste("/scratch/mjpete11/linear_models/results3/","SLC44_organ.csv",sep=""))
print(res44)
sink()

# Diagnostic plots 
m44 <- check_model(model44, show_dots=TRUE, se=FALSE)
check_normality(model44)
check_collinearity(model44)
pdf(paste("/scratch/mjpete11/linear_models/results3/","SLC44_check_model.pdf",sep=""))
m44
dev.off()

# Subset one gene from the organ df
A45 <- organs2 %>% filter(gene == "SLC25A45")
model45 <- lm(value ~ 0 + organ + SUBJID, A45)
res45 <- summary(model45)
sink(paste("/scratch/mjpete11/linear_models/results3/","SLC45_organ.csv",sep=""))
print(res45)
sink()

# Diagnostic plots 
m45 <- check_model(model45, show_dots=TRUE, se=FALSE)
check_normality(model45)
check_collinearity(model45)
pdf(paste("/scratch/mjpete11/linear_models/results3/","SLC45_check_model.pdf",sep=""))
m45
dev.off()

# Subset one gene from the organ df
A46 <- organs2 %>% filter(gene == "SLC25A46")
model46 <- lm(value ~ 0 + organ + SUBJID, A46)
res46 <- summary(model46)
sink(paste("/scratch/mjpete11/linear_models/results3/","SLC46_organ.csv",sep=""))
print(res46)
sink()

# Diagnostic plots 
m46 <- check_model(model46, show_dots=TRUE, se=FALSE)
check_normality(model46)
check_collinearity(model46)
pdf(paste("/scratch/mjpete11/linear_models/results3/","SLC46_check_model.pdf",sep=""))
m46
dev.off()

# Subset one gene from the organ df
A47 <- organs2 %>% filter(gene == "SLC25A47")
model47 <- lm(value ~ 0 + organ + SUBJID, A47)
res47 <- summary(model47)
sink(paste("/scratch/mjpete11/linear_models/results3/","SLC47_organ.csv",sep=""))
print(res47)
sink()

# Diagnostic plots 
m47 <- check_model(model47, show_dots=TRUE, se=FALSE)
check_normality(model47)
check_collinearity(model47)
pdf(paste("/scratch/mjpete11/linear_models/results3/","SLC47_check_model.pdf",sep=""))
m47
dev.off()

# Subset one gene from the organ df
A48 <- organs2 %>% filter(gene == "SLC25A48")
model48 <- lm(value ~ 0 + organ + SUBJID, A48)
res48 <- summary(model48)
sink(paste("/scratch/mjpete11/linear_models/results3/","SLC48_organ.csv",sep=""))
print(res48)
sink()

# Diagnostic plots 
m48 <- check_model(model48, show_dots=TRUE, se=FALSE)
check_normality(model48)
check_collinearity(model48)
pdf(paste("/scratch/mjpete11/linear_models/results3/","SLC48_check_model.pdf",sep=""))
m48
dev.off()

# Diagnostic plots 
m49 <- check_model(model49, show_dots=TRUE, se=FALSE)
check_normality(model49)
check_collinearity(model49)
pdf(paste("/scratch/mjpete11/linear_models/results3/","SLC49_check_model.pdf",sep=""))
m49
dev.off()

# Subset one gene from the organ df
A49 <- organs2 %>% filter(gene == "SLC25A51")
model49 <- lm(value ~ 0 + organ + SUBJID, A49)
res49 <- summary(model49)
sink(paste("/scratch/mjpete11/linear_models/results3/","SLC51_organ.csv",sep=""))
print(res49)
sink()

# Diagnostic plots 
m49 <- check_model(model49, show_dots=TRUE, se=FALSE)
check_normality(model49)
check_collinearity(model49)
pdf(paste("/scratch/mjpete11/linear_models/results3/","SLC51_check_model.pdf",sep=""))
m49
dev.off()

# Subset one gene from the organ df
A52 <- organs2 %>% filter(gene == "SLC25A52")
model52 <- lm(value ~ 0 + organ + SUBJID, A52)
res52 <- summary(model52)
sink(paste("/scratch/mjpete11/linear_models/results3/","SLC52_organ.csv",sep=""))
print(res52)
sink()

# Diagnostic plots 
m52 <- check_model(model52, show_dots=TRUE, se=FALSE)
check_normality(model52)
check_collinearity(model52)
pdf(paste("/scratch/mjpete11/linear_models/results3/","SLC52_check_model.pdf",sep=""))
m52
dev.off()

# Subset one gene from the organ df
A53 <- organs2 %>% filter(gene == "SLC25A53")
model53 <- lm(value ~ 0 + organ + SUBJID, A53)
res53 <- summary(model53)
sink(paste("/scratch/mjpete11/linear_models/results3/","SLC53_organ.csv",sep=""))
print(res53)
sink()

# Diagnostic plots 
m53 <- check_model(model53, show_dots=TRUE, se=FALSE)
check_normality(model53)
check_collinearity(model53)
pdf(paste("/scratch/mjpete11/linear_models/results3/","SLC53_check_model.pdf",sep=""))
m53
dev.off()


