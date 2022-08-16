# PCA plots of unadjusted, quantile normalize, and combat normalized GTEx count data

library(data.table)
library(stats)
library(ggplot2)

# Read in quantile normalized GTEx data
counts <- fread("/scratch/mjpete11/linear_models/data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct", sep="\t")

# PCA on unnormalized data
unnormalized_pca <- prcomp(counts[,3:ncol(counts)])
#unnormalized_pca <- prcomp(counts[,3:20])

# convert to dataframe
unnormalized_pca <- as.data.frame(unnormalized_pca[2]$rotation)

#as above, create a PCA plot for comparison to the uncorrected data
cols <- c("UHR" = "#481567FF", "HBR" = "#1F968BFF")
p1 = ggplot(data=unnormalized_pca, aes(x=PC1, y=PC2))
p1 = p1 + geom_point(size=3)
p1 = p1 + stat_ellipse(type="norm", linetype=2)
p1 = p1 + labs(title="PCA of unnormalized GTEx counts")
p1 = p1 + scale_colour_manual(values = cols)

pdf(file="/scratch/mjpete11/linear_models/results2/unnormalized_PCA.pdf")
p1
dev.off()
