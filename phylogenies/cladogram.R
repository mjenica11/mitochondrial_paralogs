library(ape)
library(ggtree)
library(ggrepel)
library(phangorn)
library(Biostrings)
library(VennDiagram)

# Import data
query_SLC25A1 <- readAAMultipleAlignment(filepath="/scratch/mjpete11/phylogenies/fastas/copy_SLC25A1_alignment.fasta",format="fasta")

query_SLC25A2 <- readAAMultipleAlignment(filepath="/scratch/mjpete11/phylogenies/fastas/copy_SLC25A2_alignment.fasta",
										 format="fasta")

query_SLC25A3 <- readAAMultipleAlignment(filepath="/scratch/mjpete11/phylogenies/fastas/copy_SLC25A3_alignment.fasta",
										 format="fasta")

# Venn diagram
venn.diagram(
   x = list(rownames(query_SLC25A1),rownames(query_SLC25A1),rownames(query_SLC25A1)),
   category.names = c("seed SLC25A1", "seed SLC25A2", "seed SLC25A3"),
   filename = "venn_diagram.png",
   output=TRUE,

   # Set names
   cat.fontface = "bold",
   cat.dist = c(0.055, 0.055, 0.055),

   # Circles
   lwd = 0.5)

# Convert alignment object to phyDat object
query_SLC25A1 <- as.phyDat(query_SLC25A1)
query_SLC25A2 <- as.phyDat(query_SLC25A2)
query_SLC25A3 <- as.phyDat(query_SLC25A3)

# Calculate distance matrix with Blosum62 matrix
dist1 <- dist.ml(query_SLC25A1, model="Blosum62")
dist2 <- dist.ml(query_SLC25A2, model="Blosum62")
dist3 <- dist.ml(query_SLC25A3, model="Blosum62")

# Make tree with Neighbor-Joining
tree1 <- NJ(dist1)
tree2 <- NJ(dist2)
tree3 <- NJ(dist3)

# Make cladogram with phangorn package
#pdf("cladogram.pdf")
#plot.phylo(tree, type="cladogram", use.edge.length=FALSE, no.margin=TRUE)
#dev.off()

# Make cladogram with ggtree package
pdf("SLC25A1_cladogram.pdf")
plt1 <- ggtree(tree1, branch.length='none', layout='circular') 
plt1 + geom_tiplab(geom="text") 
dev.off()

pdf("SLC25A2_cladogram.pdf")
plt2 <- ggtree(tree2, branch.length='none', layout='circular') 
plt2 + geom_tiplab(geom="text") 
dev.off()

pdf("SLC25A3_cladogram.pdf")
plt3 <- ggtree(tree3, branch.length='none', layout='circular') 
plt3 + geom_tiplab(geom="text") 
dev.off()
