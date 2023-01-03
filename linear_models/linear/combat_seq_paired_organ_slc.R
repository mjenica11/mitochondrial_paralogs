# Purpose: Linear models to asses if SLC paralogs are dependent on energetic state?

# Libraries
library(plyr)
#library(limma)
library(data.table)
library(tidyverse)
library(reshape)
library(stringr)
library(ggplot2)
library(ggpubr)

# Read in voom quantile normalized counts
organs <- fread("/scratch/mjpete11/linear_models/data/voom_quantile_normalized_counts3.csv", sep=",") # float

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

# Function to add multiple testing adjusted p-values for plotting

pvals <- function(GENE){
	# Subset to SLC gene
	dat <- organs %>% filter(gene == GENE)
	dat2 <- na.omit(dat)
	dat3 <- distinct(na.omit(dat2))
	# Fit the first lm; first part of two squares least regression
	model_organs <- lm(voom ~ SMRIN + SMTSISCH, data=dat3)
	# Add fitted values to dataframe
	dat3$resids <- residuals(model_organs)
	# Regress organ on SLC using previous fitted regression values as offset
	model_organs2 <- lm(voom ~ organ + resids, data=dat3)
	# Calculate adjusted p values
	obj <- summary(lm(voom ~ 0 + organ + resids, data=dat3))
	# Add p-values to dataframe
	dat3$pval <- obj$coefficients[1,4] 
	return(dat3)
}
dats <- Map(pvals, GENE=SLC)

# Set y lim
range(organs$voom) # c(0.464, 4.65) at 5 counts filtering threshold

# Density plot
density_plot <- function(GENE){
     dat <- organs %>% filter(gene == GENE)
     dat2 <- distinct(na.omit(dat))
#	 dat2$value <- round(dat2$value, 2)
     p <- ggplot(dat2, aes(voom)) +
	   geom_density() +
#	   geom_vline(xintercept=log(5), color="red") + # log of the filtering threshold
       ylab("density") +
       xlab("logCPM") +
	   coord_cartesian(xlim=c(0,5), ylim=c(0,35)) +
       ggtitle(paste0("Density plots of combat_seq voom quantile normalized ", GENE, " expression"))  
 	   ggsave(paste0("/scratch/mjpete11/linear_models/linear/voom_combat_seq_density_plots/", GENE, ".png"), p, device="png")
#     return(p)
}
plts <- Map(density_plot, GENE=SLC)

# Regress out effect of RIN and ischemic time and remove from the model 
# before making the violin plots

###### TEST ######
## Double check that for each value in the subject ID column there is a 
## corresponding heart and liver value in the organ column
#any(tapply(dA1$organ, dA1$SUBJID, function(x) length(unique(x))) != 2) #FALSE
## Are the number of rows divided by 2 equal to the number of unique subject IDs?
#nrow(dA1)/2==148 # TRUE
## Regress out the effect due to RIN and ischemic time
#model1 <- lm(value ~ SMRIN + SMTSISCH, data=dA1)
## Fit model
#dA1$resids <- residuals(model1)
## Linear regression on organ
#model2 <- lm(value ~ organ + resids, data=dA1)
## Add fitted SLC values to df for plotting
#dA1$fitted_values <- model2$fitted.values
#
## Test violin jitter plot function
#p <- ggplot(dA1, aes(x = organ, y = fitted_values)) +
#	 geom_violin(aes(colour = organ)) +
#	 geom_jitter(aes(colour = organ))  +
#	 ylab("fitted value") +
#	 xlab("organ") +
#	 ggtitle(paste0("Jitter plot of SLC25A1 expression between heart and liver after 2SRI")) +
##   scale_y_continuous(breaks = round(seq(min(dA1$fitted_values), max(dA1$fitted_values), by = 2000),0)) +
#	 stat_compare_means(method = "t.test") +
#	 stat_summary(fun.data = "mean_cl_boot", geom = "crossbar", colour = "red", width = 0.2)
#p
#ggsave(paste0("SLC25A1_violin_plot.pdf"))
####### TEST ######
#
## SLCs that passed the filtering threshold and can be plotted

library(rstatix)
# Violin and jitter plot function
violin_plots <- function(GENE, DATA){
			# Violin plot
			p <- ggplot(dat, aes(x = organ, y = voom)) +
			geom_violin(aes(colour = organ)) +
			geom_jitter(aes(colour = organ))  +
#	        coord_cartesian(ylim=c(0,25)) +
			ylim(c(0,5)) +
			ylab("logCPM") +
			xlab("organ") +
			ggtitle(paste0("Violin plot of ", GENE, " expression between heart and liver after 2SRI")) +
			stat_summary(fun.data = "mean_cl_boot", geom = "crossbar", colour = "red", width = 0.2) +
	        ggsave(paste0("/scratch/mjpete11/linear_models/linear/voom_combat_seq_violin_plots/", GENE, ".png"), p, device="png")
}
plots <- Map(violin_plots, GENE=SLC, DATA=dat)
plots[[1]]
plots
#
## Plot all of the genes next to each other
## Two squares residual inclusion on all of the SLCs
#modelA <- lm(value ~ SMRIN + SMTSISCH, data=organs)
#organs$resids <- residuals(modelA)
#modelB <- lm(value ~ organ + resids, data=organs)
#organs$fitted_values <- fitted.values(modelB)
#organs <- distinct(na.omit(organs))
#
#library(viridis)
#library(hrbrthemes)
#p2 <- organs %>% 
#		mutate(organ = factor(organ, levels = c("heart", "liver"))) %>%
#		ggplot(aes(fill=organ, y=fitted_values, x=gene)) +
#		geom_violin(position="dodge", alpha=0.5) +
#		scale_fill_viridis(discrete=T, name="") +
##		theme_ipsum() +
#		xlab("organ") +
#		ylab("counts")
#ggsave(paste0("/scratch/mjpete11/linear_models/linear/combat_seq_paired_violin_plots/SLC_grouped2.pdf"), p2, device="pdf")
