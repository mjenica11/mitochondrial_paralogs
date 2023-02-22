# Purpose: Linear models to asses if SLC paralogs are dependent on energetic state?

# Libraries
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
# Subset to the range of expected values
# I don't think there are any but my wifi sucks and I can't see the plots rn

# Remove samples >6 standard deviations away
median(organs1$fitted_values) # 2.79
mean(organs1$fitted_values) # 3.96
sd(organs1$fitted_values) # 4.00 
sd(organs1$fitted_values) * 6 # +/- 24.0 
above <- 3.96 + 24.0 # 27.96
below <- 3.96 - 24.0 # -20.04

# Are there any samples outside of this range?
range(organs1$fitted_values) # -1.62 to 29.0 

# Number of samples outside of this range
outlier_above <- organs1[organs1$fitted_values > 27.96,] # 1 sample at 29.0 log2(CPM)
outlier_below <- organs1[organs1$fitted_values < -20.04,] # 0 samples 

# Drop outliers
organs2 <- organs1[-1012,] 

# Number of samples dropped should be equal to 97
nrow(organs1) - nrow(organs2)== 1 # TRUE

rm(violin_plots)
# Violin and jitter plot function
violin_plots <- function(GENE){
	dat <- organs1 %>% filter(gene == GENE)
    dat2 <- na.omit(dat)
    dat3 <- distinct(na.omit(dat2))
    # Regress out effect of RIN and ischemic time on SLC
    model1 <- lm(voom ~ SMRIN + SMTSISCH, data=dat3)
	# Add fitted values to dataframe
	dat3$resids <- residuals(model1)
	# Regress organ on SLC using previous fitted regression values as offset
	model2 <- lm(voom ~ organ + resids, data=dat3)
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
	#annotate(geom = "text", x = 1.5, y = 19, label=paste0("p value: ",corrected_pval)) +
	annotate(geom = "text", x = 1.5, y = 34, label=paste0("p value: ",corrected_pval)) +
	#ylim(c(0,30)) +
	ylab("log2(CPM)") +
	xlab("organ") +
	scale_y_continuous(limits = c(-35, 35), expand = c(0,0), breaks = seq(-35, 35, by = 1)) +
#	scale_y_continuous(breaks=scales::pretty_breaks(n=20)) + # Add more y-axis tick marks
	ggtitle(paste0("Violin plot of ", GENE, " expression between heart and liver after 2SRI")) 
	ggsave(paste0("/scratch/mjpete11/linear_models/linear/voom_combat_seq_violin_plots/", GENE, ".png"), p, device="png")
}
plots <- Map(violin_plots, GENE=SLC)
violin_plots(GENE="SLC25A1", DATA=organs)
plots[[1]]
plots
