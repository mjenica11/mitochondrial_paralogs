# Purpose: Linear models to asses if SLC paralogs are dependent on energetic state?

# Libraries
library(pscl)
library(ggtext)
library(MASS)
library(plyr)
library(dplyr)
library(data.table)
library(tidyverse)
library(reshape)
library(stringr)
library(ggplot2)
library(ggpubr)
library(rstatix)
library(stats)
library(scales)

# Read in voom quantile normalized counts
organs <- fread("/scratch/mjpete11/linear_models/data/voom_qnorm_counts1.csv", sep=",") # float

# Drop the index column
organs$V1 <- NULL

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

# Set y lim
range(organs$voom) # c(0.464, 31.73) at 5 counts filtering threshold

# Density plot
#density_plot <- function(GENE){
#     dat <- organs %>% filter(gene == GENE)
#     dat2 <- distinct(na.omit(dat))
#     p <- ggplot(dat2, aes(voom)) +
#	 geom_density() +
#     ylab("density") +
#     xlab("logCPM") +
#	 coord_cartesian(xlim=c(0,35), ylim=c(0,20)) +
#     ggtitle(paste0("Density plots of combat_seq voom quantile normalized ", GENE, " expression"))  
# 	 ggsave(paste0("/scratch/mjpete11/linear_models/linear/voom_qnorm_density_plots3/", GENE, ".png"), p, device="png")
#}
#plts <- Map(density_plot, GENE=SLC)

# Regress out effect of RIN and ischemic time and remove from the model 
# before making the violin plots

##### TEST #####
contained <- list()
for (paralog in SLC){
	datt <- organs %>% filter(gene == paralog)
	datt2 <- na.omit(datt)
	datt3 <- distinct(na.omit(datt2))
	# Regress out effect of RIN and ischemic time on SLC
	model1 <- lm(voom ~ SMRIN + SMTSISCH, data=datt3)
	# Add fitted values to dattaframe
	datt3$resids <- residuals(model1)
	# Regress organ on SLC using previous fitted regression values as offset
	model2 <- lm(voom ~ organ + resids, data=datt3)
	# Add fitted values from second linear model to dattaframs
	datt3$fitted_values <- fitted.values(model2)
	contained[[paralog]] <- datt3
}
organs1 <- ldply(contained, data.frame)  


##### TEST #####

# Write the striated samples to file
# There are none

# Remove samples >6 standard deviations away
median(organs1$counts) # 801
mean(organs1$counts) # 3356
sd(organs1$counts) # 9134 
sd(organs1$counts) * 6 # +/- 54804 
above <- 3356 + 54804 # 58160
below <- 3356 - 54804 # -51448 but we will cap at 0

# Are there any samples outside of this range?
range(organs1$counts) # 0 to 229,728 

# Number of samples outside of this range
outlier_above <- organs1[organs1$counts > 58160,] # 100 samples 
outlier_below <- organs1[organs1$counts < 0,] # 0 samples 

# Drop outliers
organs2 <- organs1[!(row.names(organs1) %in% row.names(outlier_above)),]

# Number of samples dropped should be equal to 97
nrow(organs1) - nrow(organs2)== 100 # TRUE

# How many zeros were present before applying voom() but after filtering?
nrow(organs1[organs1$counts==0,]) # 1262

# Add an additional column where voom values are the same except
# for genes that were originally zero, convert back to zero
# voom() add a value of 2 for zero values to avoid taking the log of zero
organs2$voom[organs2$counts==0] <- 0

# How many zeros are present in the voom transformed values column? 
nrow(organs2[organs2$voom==0,]) # 1262

rm(ZINB_violin_plots)
# Zero-inflated negative binomial 2SRI violin and jitter plot function
ZINB_violin_plots <- function(GENE){
	dat <- organs2 %>% filter(gene == GENE)
    dat2 <- na.omit(dat)
    dat3 <- distinct(na.omit(dat2))
    # Regress out effect of RIN and ischemic time on SLC
    model1 <- zeroinfl(counts ~ SMRIN + SMTSISCH, data=dat3)
	# Add fitted values to dataframe
	dat3$resids <- residuals(model1)
	# Regress organ on SLC using previous fitted regression values as offset
	model2 <- zeroinfl(counts ~ organ + resids, data=dat3)
	# Add fitted values from second linear model to dataframs
	dat3$fitted_values <- fitted.values(model2)
	p_val <- wilcox.test(formula=fitted_values~organ,data=dat3)$p.value
	corrected_pval <- p.adjust(p_val, method="bonferroni", n=53)
    # Violin plot
	p <- ggplot(dat3, aes(x = organ, y = fitted_values, fill=organ)) +
	geom_violin(width=1, position=position_dodge(width=0.5)) +
	scale_fill_manual(values=c("lightgreen", "purple")) +
	geom_jitter(width=0.3) +
	stat_summary(fun.data="mean_sdl", geom="crossbar", width=0.1, alpha=0.1) +
	annotate(geom = "text", x = 1.5, y = max(dat3$counts)+5, label=paste0("p value: ",corrected_pval)) +
	#ylim(c(0,30)) +
	ylab("counts") +
	xlab("organ") +
#	scale_y_continuous(limits = c(-35, 35), expand = c(0,0), breaks = seq(-35, 35, by = 1)) +
	theme(plot.title = element_textbox_simple(size=10)) + # automatically wrap long titles
	ggtitle(paste0("Violin plot of ", GENE, " expression between heart and liver after ZINB 2SRI")) 
	# save to different directories depending on if y ticks or not
#	ggsave(paste0("/scratch/mjpete11/linear_models/linear/ZINB_2SRI_violin_plots/no_yticks/", GENE, ".png"), p, device="png")
	ggsave(paste0("/scratch/mjpete11/linear_models/linear/ZINB_2SRI_violin_plots/yticks/", GENE, ".png"), p, device="png")
}
# Have to plot each SLC25 individually since some are zero inflated and some aren't...
ZINB_violin_plots(GENE="SLC25A2")
ZINB_violin_plots(GENE="SLC25A6")
ZINB_violin_plots(GENE="UCP1")
ZINB_violin_plots(GENE="SLC25A21")
ZINB_violin_plots(GENE="SLC25A41")
ZINB_violin_plots(GENE="SLC25A47")
ZINB_violin_plots(GENE="SLC25A48")
ZINB_violin_plots(GENE="SLC25A52")

# NB violin plots
# MASS is loaded by pscl, so I have to unload pscl to use MASS::glm.nb
unloadNamespace("pscl")

rm(NB_violin_plots)
# Negative binomial 2SRI violin and jitter plot function
NB_violin_plots <- function(GENE){
	dat <- organs2 %>% filter(gene == GENE)
    dat2 <- na.omit(dat)
    dat3 <- distinct(na.omit(dat2))
    # Regress out effect of RIN and ischemic time on SLC
    model1 <- glm.nb(counts ~ SMRIN + SMTSISCH, data=dat3)
	# Add fitted values to dataframe
	dat3$resids <- residuals(model1)
	# Regress organ on SLC using previous fitted regression values as offset
	model2 <- glm.nb(counts ~ organ + resids, data=dat3)
	# Add fitted values from second linear model to dataframs
	dat3$fitted_values <- fitted.values(model2)
	p_val <- wilcox.test(formula=fitted_values~organ,data=dat3)$p.value
	corrected_pval <- p.adjust(p_val, method="bonferroni", n=53)
    # Violin plot
	p <- ggplot(dat3, aes(x = organ, y = fitted_values, fill=organ)) +
	geom_violin(width=1, position=position_dodge(width=0.5)) +
	scale_fill_manual(values=c("lightgreen", "purple")) +
	geom_jitter(width=0.3) +
	stat_summary(fun.data="mean_sdl", geom="crossbar", width=0.1, alpha=0.1) +
	annotate(geom = "text", x = 1.5, y = max(dat3$counts)+5, label=paste0("p value: ",corrected_pval)) +
	#ylim(c(0,30)) +
	ylab("counts") +
	xlab("organ") +
#	scale_y_continuous(limits = c(-35, 35), expand = c(0,0), breaks = seq(-35, 35, by = 1)) +
	theme(plot.title = element_textbox_simple(size=10)) + # automatically wrap long titles
	ggtitle(paste0("Violin plot of ", GENE, " expression between heart and liver after NB 2SRI")) 
	ggsave(paste0("/scratch/mjpete11/linear_models/linear/ZINB_2SRI_violin_plots/no_yticks/", GENE, ".png"), p, device="png")
}
NB_violin_plots(GENE="SLC25A1")
NB_violin_plots(GENE="SLC25A3")
NB_violin_plots(GENE="SLC25A4")
NB_violin_plots(GENE="SLC25A5")
NB_violin_plots(GENE="UCP2")
NB_violin_plots(GENE="UCP3")
NB_violin_plots(GENE="SLC25A10")
NB_violin_plots(GENE="SLC25A11")
NB_violin_plots(GENE="SLC25A12")
NB_violin_plots(GENE="SLC25A13")
NB_violin_plots(GENE="SLC25A14")
NB_violin_plots(GENE="SLC25A15")
NB_violin_plots(GENE="SLC25A16")
NB_violin_plots(GENE="SLC25A17")
NB_violin_plots(GENE="SLC25A18")
NB_violin_plots(GENE="SLC25A19")
NB_violin_plots(GENE="SLC25A20")
NB_violin_plots(GENE="SLC25A22")
NB_violin_plots(GENE="SLC25A23")
NB_violin_plots(GENE="SLC25A24")
NB_violin_plots(GENE="SLC25A25")
NB_violin_plots(GENE="SLC25A26")
NB_violin_plots(GENE="SLC25A27")
NB_violin_plots(GENE="SLC25A28")
NB_violin_plots(GENE="SLC25A29")
NB_violin_plots(GENE="SLC25A30")
NB_violin_plots(GENE="SLC25A31")
NB_violin_plots(GENE="SLC25A32")
NB_violin_plots(GENE="SLC25A33")
NB_violin_plots(GENE="SLC25A34")
NB_violin_plots(GENE="SLC25A35")
NB_violin_plots(GENE="SLC25A36")
NB_violin_plots(GENE="SLC25A37")
NB_violin_plots(GENE="SLC25A38")
NB_violin_plots(GENE="SLC25A39")
NB_violin_plots(GENE="SLC25A40")
NB_violin_plots(GENE="SLC25A42")
NB_violin_plots(GENE="SLC25A43")
NB_violin_plots(GENE="SLC25A44")
NB_violin_plots(GENE="SLC25A45")
NB_violin_plots(GENE="SLC25A46")
NB_violin_plots(GENE="MTCH1")
NB_violin_plots(GENE="MTCH2")
NB_violin_plots(GENE="SLC25A52")
NB_violin_plots(GENE="SLC25A53")

