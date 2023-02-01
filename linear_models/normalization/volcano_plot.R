# Add the estimated gene expression values to the organs metadataframe

library(RVA)
library(data.table)
library(tidyverse)
library(ggplot2)

# Read in topTable results from the batch_voom_qnorm.R script
tab <- fread("/scratch/mjpete11/linear_models/data/batch_voom_qnorm_topTable.csv")

# Add the posterior gene expression estimates to the counts df
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

# In limma, I subsetted the eBayes() results with topTable() by setting
# the genelist parameter to "SLC", but the output only contains the row
# index and not the gene name
# The row index of the count matrix supplied as input to voom() is in the
# same order as the organs metadata frame, so I will use index and gene names column from the filtered counts matrix to assign a gene name to each 
# row index in the topTable() results
counts <- fread("/scratch/mjpete11/linear_models/data/filtered_counts.csv", sep=",", select=c(1,2))
nrow(counts)==nrow(tab) #FALSE
tab2 <- merge(counts, tab, by="V1")  
nrow(counts)==nrow(tab2) #FALSE

# Range of logFC and p values
range(tab$logFC) # -14.32 to 17.21
range(tab$P.Value) # 0.00 to 0.048 
length(which(tab$P.Value==0))

# rename the moderated t statistic column so the interpreter does not confused
# it with the t() function
colnames(tab2)[7] <- c("moderated_t")
head(tab2)

# Subset just the rows with the SLC values to check they are all there
tab3 <- tab2[tab2$Description %in% SLC,]

# Drop the index column in the topTable() results
tab$V1 <- NULL

# Convert the gene name column the rownames
rownames(tab) <- make.names(tab2$Description, unique=TRUE)

# Drop the confidence intervals
tab$CI.L <- NULL 
tab$CI.R <- NULL 

# Generate volcano plot
plt <- plot_volcano(
				data = tab,
				comp.names = NULL,
				geneset = NULL,
			    geneset.FCflag = "logFC",
			    highlight.1 = SLC,
				highlight.2 = SLC,
			    upcolor = "#FF0000",
			    downcolor = "#0000FF",
				plot.save.to = "/scratch/mjpete11/linear_models/normalization/batch_voom_violin_plots/volcano.pdf",
				xlim = c(-15,18),
				ylim = c(0,1000),
			    FCflag = "logFC",
				FDRflag = "adj.P.Val",
				highlight.FC.cutoff = 1.5,
				highlight.FDR.cutoff = 0.05,
				title = "Volcano plot after blocking batch with voom()",
				xlab = "log2 Fold Change",
				ylab = "log10(FDR)")
plt								 
