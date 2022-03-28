# Convert percent identity matrix text file to a triangle 
# Read in PIM file
pim <- read.csv("pim.txt", header=FALSE, sep="", skip=1)
pim <- as.matrix(pim[,-1])

# Drop top half of matrix
pim[upper.tri(pim, diag=TRUE)] <- 0
pim <- as.matrix(pim)
dim(pim)
head(pim)

# Change the column and row names to the HGNC names
HGNC <- c("SLC25A22", "SLC25A44", "SLC25A33", "SLC25A26", "SLC25A4", "SLC25A38",
		  "SLC25A30", "SLC25A46", "SLC25A24", "SLC25A18", "SLC25A12", "SLC25A13", 
		  "SLC25A43", "SLC25A40", "SLC25A14", "SLC25A8", "SLC25A27", "SLC25A21", 
		  "SLC25A3", "SLC25A36", "SLC25A32", "SLC25A7", "SLC25A50", "SLC25A20", 
		  "SLC25A5", "SLC25A17", "SLC25A16", "SLC25A31", "SLC25A53", "SLC25A51", 
		  "SLC25A52", "SLC25A23", "SLC25A34", "SLC25A48", "SLC25A25", "SLC25A42", 
		  "SLC25A10", "SLC25A49", "SLC25A29", "SLC25A9", "SLC25A47", "SLC25A6", 
		  "SLC25A41", "SLC25A11", "SLC25A1", "SLC25A35", "SLC25A45", "SLC25A15",
		  "SLC25A2", "SLC25A39", "SLC25A28", "SLC25A19", "SLC25A37")

row.names(pim) <- HGNC 
colnames(pim) <- HGNC

# write distance matrix
write.csv(pim, "triangle_pim.txt")

# make heatmap plot
pdf("heatmap.pdf")
heatmap(pim)
dev.off()
