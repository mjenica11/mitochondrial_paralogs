# Convert percent identity matrix text file to a triangle 
# Read in PIM file
PGAM_pim <- read.csv("PGAM_pim.txt", header=FALSE, sep="", skip=1)
PGAM_pim <- as.matrix(PGAM_pim[,-1])

# Drop top half of matrix
#pim[upper.tri(pim, diag=TRUE)] <- 0
PGAM_pim <- as.matrix(PGAM_pim)
dim(PGAM_pim)
head(PGAM_pim)

PGAM_genes <- c("PGAM", "GMCL2", "c10orf120", "BTBD16")
row.names(PGAM_pim) <- PGAM_genes 
colnames(PGAM_pim) <- PGAM_genes

# make heatmap plot
pdf("PGAM_heatmap.pdf")
heatmap(PGAM_pim, margins=c(10,10))
dev.off()
