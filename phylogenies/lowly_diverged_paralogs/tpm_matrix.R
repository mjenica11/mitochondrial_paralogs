# Convert percent identity matrix text file to a triangle 
# Read in PIM file
GMCL1_pim <- read.csv("GMCL1_pim.txt", header=FALSE, sep="", skip=1)
GMCL1_pim <- as.matrix(GMCL1_pim[,-1])

# Drop top half of matrix
#pim[upper.tri(pim, diag=TRUE)] <- 0
GMCL1_pim <- as.matrix(GMCL1_pim)
dim(GMCL1_pim)
head(GMCL1_pim)

GMCL1_genes <- c("GMCL1", "GMCL2", "c10orf120", "BTBD16")
row.names(GMCL1_pim) <- GMCL1_genes 
colnames(GMCL1_pim) <- GMCL1_genes

# make heatmap plot
pdf("GMCL1_heatmap.pdf")
heatmap(GMCL1_pim, margins=c(10,10))
dev.off()
