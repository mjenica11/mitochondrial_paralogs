# Purpose: Linear models to asses if SLC paralogs are dependent on energetic state?

# Libraries
library(pscl)
library(MASS)
library(plyr)
library(dplyr)
library(data.table)
library(tidyverse)
library(reshape)
library(stringr)
library(ggplot2)
library(ggpubr)
library(ggtext)
library(rstatix)
library(stats)
library(scales)
library(edgeR)

# Read in combat_seq normalized counts 
counts <- fread("/scratch/mjpete11/linear_models/data/combat_seq_filtered.csv", sep=",")

# Drop the index column
counts$V1 <- NULL

# Generate the design matrix
# Samples are rows and columns are covariates of interest (heart and liver)

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

# Subset the SLC25 genes
counts <- as.data.frame(counts)
sub_df <- counts[counts$"Description" %in% SLC, ]

# Write to file
#write.table(sub_df, "/scratch/mjpete11/linear_models/data/SLC_df_voom_combat_seq.csv", sep=",")

# Read in count df of just SLC genes
sub_df <- read.csv("/scratch/mjpete11/linear_models/data/SLC_df_voom_combat_seq.csv", sep=",")
names(sub_df) <- gsub(x=names(sub_df), pattern="\\.", replacement="-")

# All present!
setdiff(SLC, sub_df$'Description') 

# Read in sample attributes files
file <- read.csv("/scratch/mjpete11/linear_models/data/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt", sep="\t")

# Subset heart left ventricle and liver samples into two separate dfs
heart <- file[file$SMTSD %in% "Heart - Left Ventricle", ]
liver <- file[file$SMTSD %in% "Liver", ]

# Subset the SLC25 gene count df by the sample IDs that match the IDs in the 
# GTEx sample annotation df
heart2 <- sub_df %>% select(contains(heart$SAMPID))
liver2 <- sub_df %>% select(contains(liver$SAMPID))

# Append the gene name
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

# Any duplicated rows?
any(duplicated(heart4)==TRUE) # FALSE
any(duplicated(liver4)==TRUE) # FALSE

# Change the name of the 'variable' column to 'SAMPID' to match columns
colnames(heart4)[3] <- "SAMPID" 
colnames(liver4)[3] <- "SAMPID" 

# Add a subject ID column
heart4$SUBJID <- str_sub(heart4$SAMPID, start=1L, end=10L)
liver4$SUBJID <- str_sub(liver4$SAMPID, start=1L, end=10L)

# Add column with the expression batch ID (SMGEBTCH) and the type of genotype
# or expression batch (SMGEBTCHT), the RIN number (SMRIN) and tissue type (SMTSISCH)
heart4$SMGEBTCH <- file$SMGEBTCH[match(heart4$SAMPID, file$SAMPID)]
liver4$SMGEBTCH <- file$SMGEBTCH[match(liver4$SAMPID, file$SAMPID)]

heart4$SMGEBTCHT <- file$SMGEBTCHT[match(heart4$SAMPID, file$SAMPID)]
liver4$SMGEBTCHT <- file$SMGEBTCHT[match(liver4$SAMPID, file$SAMPID)]

heart4$SMRIN <- file$SMRIN[match(heart4$SAMPID, file$SAMPID)]
liver4$SMRIN <- file$SMRIN[match(liver4$SAMPID, file$SAMPID)]

heart4$SMTSISCH <- file$SMTSISCH[match(heart4$SAMPID, file$SAMPID)]
liver4$SMTSISCH <- file$SMTSISCH[match(liver4$SAMPID, file$SAMPID)]

# Dataframe with only samples from people who donated both a heart and liver: organs
# e.g. paired samples
organs_tmp <- merge(heart4, liver4, by="SUBJID") 
length(unique(organs_tmp$SUBJID)) # 148 individuals
tmp1 <- organs_tmp[,c("SUBJID", "gene.x", "organ.x", "SAMPID.x", "value.x", "SMRIN.x", "SMTSISCH.x", "SMGEBTCH.x", "SMGEBTCHT.x")]
tmp2 <- organs_tmp[,c("SUBJID", "gene.y", "organ.y", "SAMPID.y", "value.y", "SMRIN.y", "SMTSISCH.y", "SMGEBTCH.x", "SMGEBTCHT.x")]
colnames(tmp1) <- c("SUBJID", "gene", "organ", "SAMPID", "value", "SMRIN", "SMTSISCH", "SMGEBTCH", "SMGEBTCHT")
colnames(tmp2) <- c("SUBJID", "gene", "organ", "SAMPID", "value", "SMRIN", "SMTSISCH", "SMGEBTCH", "SMGEBTCHT")
organs0 <- rbind(tmp1, tmp2)

# Are any rows duplicated?
duplicated(organs0) # Yes

# Drop duplicated rows
organs <- organs0[!duplicated(organs0),]
length(unique(organs_tmp$SUBJID)) # 148 individuals

# Change the name of 'value' column to 'combat_seq_counts'
colnames(organs)[5] <- "combat_seq_counts"

# Convert combat_seq adjusted counts to log2 and CPM
# Values of zero will be set to 2
organs$log2_cpm <- cpm(organs$combat_seq_counts, log=TRUE, prior.count=2)

# Write organs df to file
#write.table(organs, "/scratch/mjpete11/linear_models/data/organs_combat_seq.csv", sep=",")

# Drop samples from the counts df that are not from paired heart or liver
keeps <- as.vector(organs$SAMPID)
counts2 <- subset(counts, select=keeps) 

# Are there the same number of samples in the counts df and the organs metadata df? 
ncol(counts2)==nrow(organs) # TRUE

# Set y lim
range(organs$combat_seq_counts) # c(0 229728) 
range(organs$log2_cpm) # c(-4.72, 12.09) 

# Regress out effect of RIN and ischemic time and remove from the model 
# before making the violin plots

##### TEST #####
contained <- list()
for (paralog in SLC){
	datt <- organs %>% filter(gene == paralog)
	datt2 <- na.omit(datt)
	datt3 <- distinct(na.omit(datt2))
	# Regress out effect of RIN and ischemic time on SLC
	model1 <- lm(log2_cpm ~ SMRIN + SMTSISCH, data=datt3)
	# Add fitted values to dattaframe
	datt3$resids <- residuals(model1)
	# Regress organ on SLC using previous fitted regression values as offset
	model2 <- lm(log2_cpm ~ organ + resids, data=datt3)
	# Add fitted values from second linear model to dattaframs
	datt3$fitted_values <- fitted.values(model2)
	contained[[paralog]] <- datt3
}
organs1 <- ldply(contained, data.frame)  
head(organs1)
##### TEST #####

# Write the striated samples to file
# There are none

# Remove samples >6 standard deviations away
median(organs1$fitted_values) # 3.94
mean(organs1$fitted_values) # 3.29
sd(organs1$fitted_values) # 3.52 
sd(organs1$fitted_values) * 6 # +/- 21.12 
above <- 3.29 + 21.12 # 24.41 
below <- 3.29 - 21.12 # -17.83

# Are there any samples outside of this range?
range(organs1$fitted_values) # -5.56, 11.9 --> No 

# NB violin plots
# MASS is loaded by pscl, so I have to unload pscl to use MASS::glm.nb
unloadNamespace("pscl")

# Poisson 2SRI model
rm(plots)
rm(gaussian_2SRI_violin_plots)
# Negative binomial 2SRI violin and jitter plot function
gaussian_2SRI_violin_plots <- function(GENE){
	dat <- organs1 %>% filter(gene == GENE)
    dat2 <- na.omit(dat)
    dat3 <- distinct(na.omit(dat2))
    # Regress out effect of RIN and ischemic time on SLC
    model1 <- glm(log2_cpm ~ SMRIN + SMTSISCH, data=dat3)
	# Add fitted values to dataframe
	dat3$resids <- residuals(model1)
	# Regress organ on SLC using previous fitted regression values as offset
	model2 <- glm(log2_cpm ~ organ + resids, data=dat3)
	# Add fitted values from second linear model to dataframs
	dat3$fitted_values <- fitted.values(model2)
#	p_val <- wilcox.test(formula=fitted_values~organ,data=dat3)$p.value
#	corrected_pval <- p.adjust(p_val, method="bonferroni", n=53)
    # Violin plot
	p <- ggplot(dat3, aes(x = organ, y = fitted_values, fill=organ)) +
				stat_compare_means(method = "wilcox.test", 
								   aes(label = paste("adj.p_value =", after_stat(!!str2lang("p.adj"))*53)), 
								   label.x = 1.25, 
								   label.y = max(dat3[["fitted_values"]]) + 1.5,
								   paired = TRUE) +
				geom_violin(width=1, position=position_dodge(width=0.5)) +
				scale_fill_manual(values=c("lightgreen", "purple")) +
				geom_jitter(width=0.3) +
				stat_summary(fun.data="mean_sdl", geom="crossbar", width=0.1, alpha=0.1) +
				ylab("log2_cpm") +
				xlab("organ") +
				theme(plot.title = element_textbox_simple(size=10)) + # automatically wrap long titles
				ggtitle(paste0("Violin plot of ", GENE, " expression between heart and liver after 2SRI")) 
				ggsave(paste0("/scratch/mjpete11/linear_models/linear/2SRI/gaussian_2SRI_violin_plots/", GENE, ".png"), p, device="png")
}
plots <- Map(gaussian_2SRI_violin_plots, GENE=SLC)



##################### OLD ###################################################
#rm(plots)
#rm(NB_violin_plots)
# Negative binomial 2SRI violin and jitter plot function
#NB_violin_plots <- function(GENE){
#	dat <- organs2 %>% filter(gene == GENE)
#    dat2 <- na.omit(dat)
#    dat3 <- distinct(na.omit(dat2))
    # Regress out effect of RIN and ischemic time on SLC
#    model1 <- glm.nb(combat_seq_counts ~ SMRIN + SMTSISCH, data=dat3)
	# Add fitted values to dataframe
#	dat3$resids <- residuals(model1)
	# Regress organ on SLC using previous fitted regression values as offset
#	model2 <- glm.nb(combat_seq_counts ~ organ + resids, data=dat3)
	# Add fitted values from second linear model to dataframs
#	dat3$fitted_values <- fitted.values(model2)
#	p_val <- wilcox.test(formula=fitted_values~organ,data=dat3)$p.value
#	corrected_pval <- p.adjust(p_val, method="bonferroni", n=53)
    # Violin plot
#	p <- ggplot(dat3, aes(x = organ, y = fitted_values, fill=organ)) +
#	geom_violin(width=1, position=position_dodge(width=0.5)) +
#	scale_fill_manual(values=c("lightgreen", "purple")) +
#	geom_jitter(width=0.3) +
#	stat_summary(fun.data="mean_sdl", geom="crossbar", width=0.1, alpha=0.1) +
#	annotate(geom = "text", x = 1.5, y = 34, label=paste0("p value: ",corrected_pval)) +
	#ylim(c(0,30)) +
#	ylab("combat_seq_counts") +
#	xlab("organ") +
#	scale_y_continuous(limits = c(0, 73898), expand = c(0,0), breaks = seq(0, 73898, by = 50)) +
#	theme(plot.title = element_textbox_simple(size=10)) + # automatically wrap long titles
#	ggtitle(paste0("Violin plot of ", GENE, " expression between heart and liver after NB 2SRI")) 
#	ggsave(paste0("/scratch/mjpete11/linear_models/linear/ZINB_2SRI_violin_plots/", GENE, ".png"), p, device="png")
#}
#plots <- Map(NB_violin_plots, GENE=SLC)

#rm(ZINB_violin_plots)
# Violin and jitter plot function
#ZINB_violin_plots <- function(GENE){
#	dat <- organs2 %>% filter(gene == GENE)
    #dat2 <- na.omit(dat)
    # Regress out effect of RIN and ischemic time on SLC
 #   model1 <- zeroinfl(combat_seq_counts ~ SMRIN + SMTSISCH, data=dat)
	# Add fitted values to dataframe
#	dat$resids <- residuals(model1)
	# Regress organ on SLC using previous fitted regression values as offset
#	model2 <- zeroinfl(combat_seq_counts ~ organ + resids, data=dat)
	# Add fitted values from second linear model to dataframs
#	dat$fitted_values <- fitted.values(model2)
#	p_val <- wilcox.test(formula=fitted_values~organ,data=dat)$p.value
#	corrected_pval <- p.adjust(p_val, method="bonferroni", n=53)
    # Violin plot
#	p <- ggplot(dat, aes(x = organ, y = fitted_values, fill=organ)) + geom_violin(width=1, position=position_dodge(width=0.5)) +
#	scale_fill_manual(values=c("lightgreen", "purple")) +
#	geom_jitter(width=0.3) +
#	stat_summary(fun.data="mean_sdl", geom="crossbar", width=0.1, alpha=0.1) +
	#annotate(geom = "text", x = 1.5, y = 19, label=paste0("p value: ",corrected_pval)) +
#	annotate(geom = "text", x = 1.5, y = 34, label=paste0("p value: ",corrected_pval)) +
	#ylim(c(0,30)) +
#	ylab("log2(CPM)") +
#	xlab("organ") +
#	scale_y_continuous(limits = c(-35, 35), expand = c(0,0), breaks = seq(-35, 35, by = 1)) +
#	scale_y_continuous(breaks=scales::pretty_breaks(n=20)) + # Add more y-axis tick marks
#	ggtitle(paste0("Violin plot of ", GENE, "expression between heart and liver after ZINB 2SRI")) 
#	ggsave(paste0("/scratch/mjpete11/linear_models/linear/ZINB_2SRI_violin_plots/", GENE, ".png"), p, device="png")
#}
#plots <- Map(violin_plots, GENE=SLC)

# Have to plot each SLC25 individually since some are zero inflated and some aren't...
#ZINB_violin_plots(GENE="SLC25A2")
#ZINB_violin_plots(GENE="SLC25A6")
#ZINB_violin_plots(GENE="UCP1")
#ZINB_violin_plots(GENE="SLC25A21")
#ZINB_violin_plots(GENE="SLC25A41")
#ZINB_violin_plots(GENE="SLC25A47")
#ZINB_violin_plots(GENE="SLC25A48")
#ZINB_violin_plots(GENE="SLC25A52")


#NB_violin_plots(GENE="SLC25A1")
#NB_violin_plots(GENE="SLC25A3")
#NB_violin_plots(GENE="SLC25A4")
#NB_violin_plots(GENE="SLC25A5")
#NB_violin_plots(GENE="UCP2")
#NB_violin_plots(GENE="UCP3")
#NB_violin_plots(GENE="SLC25A10")
#NB_violin_plots(GENE="SLC25A11")
#NB_violin_plots(GENE="SLC25A12")
#NB_violin_plots(GENE="SLC25A13")
#NB_violin_plots(GENE="SLC25A14")
#NB_violin_plots(GENE="SLC25A15")
#NB_violin_plots(GENE="SLC25A16")
#NB_violin_plots(GENE="SLC25A17")
#NB_violin_plots(GENE="SLC25A18")
#NB_violin_plots(GENE="SLC25A19")
#NB_violin_plots(GENE="SLC25A20")
#NB_violin_plots(GENE="SLC25A22")
#NB_violin_plots(GENE="SLC25A23")
#NB_violin_plots(GENE="SLC25A24")
#NB_violin_plots(GENE="SLC25A25")
#NB_violin_plots(GENE="SLC25A26")
#NB_violin_plots(GENE="SLC25A27")
#NB_violin_plots(GENE="SLC25A28")
#NB_violin_plots(GENE="SLC25A29")
#NB_violin_plots(GENE="SLC25A30")
#NB_violin_plots(GENE="SLC25A31")
#NB_violin_plots(GENE="SLC25A32")
#NB_violin_plots(GENE="SLC25A33")
#NB_violin_plots(GENE="SLC25A34")
#NB_violin_plots(GENE="SLC25A35")
#NB_violin_plots(GENE="SLC25A36")
#NB_violin_plots(GENE="SLC25A37")
#NB_violin_plots(GENE="SLC25A38")
#NB_violin_plots(GENE="SLC25A39")
#NB_violin_plots(GENE="SLC25A40")
#NB_violin_plots(GENE="SLC25A42")
#NB_violin_plots(GENE="SLC25A43")
#NB_violin_plots(GENE="SLC25A44")
#NB_violin_plots(GENE="SLC25A45")
#NB_violin_plots(GENE="SLC25A46")
#NB_violin_plots(GENE="MTCH1")
#NB_violin_plots(GENE="MTCH2")
#NB_violin_plots(GENE="SLC25A52")
#NB_violin_plots(GENE="SLC25A53")

