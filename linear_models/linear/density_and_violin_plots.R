# Purpose: Linear models to asses if SLC paralogs are dependent on energetic state?

# Libraries
library(plyr)
library(data.table)
library(tidyverse)
library(reshape)
library(stringr)
library(ggplot2)
library(ggpubr)
library(rstatix)

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
density_plot <- function(GENE){
     dat <- organs %>% filter(gene == GENE)
     dat2 <- distinct(na.omit(dat))
     p <- ggplot(dat2, aes(voom)) +
	 geom_density() +
     ylab("density") +
     xlab("logCPM") +
	 coord_cartesian(xlim=c(0,35), ylim=c(0,20)) +
     ggtitle(paste0("Density plots of combat_seq voom quantile normalized ", GENE, " expression"))  
 	 ggsave(paste0("/scratch/mjpete11/linear_models/linear/voom_qnorm_density_plots/", GENE, ".png"), p, device="png")
}
plts <- Map(density_plot, GENE=SLC)

# Regress out effect of RIN and ischemic time and remove from the model 
# before making the violin plots

# Violin and jitter plot function
violin_plots <- function(GENE){
	dat <- organs %>% filter(gene == GENE)
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
	P_VAL=wilcox.test(formula=fitted_values~organ,data=dat3)$p.value*53
    # Violin plot
	p <- ggplot(dat3, aes(x = organ, y = fitted_values)) +
	geom_violin(aes(colour = organ)) +
	geom_jitter(aes(colour = organ))  +
	annotate(geom = "text",x = 1.5,y = max(dat3$fitted_values),label=paste0("p value: ",P_VAL))+
	#ylim(c(0,30)) +
	ylab("logCPM") +
	xlab("organ") +
	ggtitle(paste0("Violin plot of ", GENE, " expression between heart and liver after 2SRI")) +
    stat_summary(fun.data = "mean_cl_boot", geom = "crossbar", colour = "red", width = 0.2)
	ggsave(paste0("/scratch/mjpete11/linear_models/linear/voom_qnorm_violin_plots2/", GENE, ".pdf"), p, device="pdf")
}
plots <- Map(violin_plots, GENE=SLC)
violin_plots(GENE="SLC25A1", DATA=organs)
plots[[1]]
plots
