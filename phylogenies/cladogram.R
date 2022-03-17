library(ape)
library(ggtree)
library(ggrepel)
library(phangorn)
library(Biostrings)
library(VennDiagram)

# Import data
#query_SLC25A1 <- readAAMultipleAlignment(filepath="/scratch/mjpete11/phylogenies/fastas/copy_SLC25A1_alignment.fasta",format="fasta")

query_SLC25A2 <- readAAMultipleAlignment(filepath="/scratch/mjpete11/phylogenies/fastas/copy_SLC25A2_alignment.fasta",
										 format="fasta")

#query_SLC25A3 <- readAAMultipleAlignment(filepath="/scratch/mjpete11/phylogenies/fastas/copy_SLC25A3_alignment.fasta",
#										 format="fasta")

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
query_SLC25A2 <- as.phyDat(query_SLC25A2)

# Calculate distance matrix with Blosum62 matrix
dist_obj <- dist.ml(query_SLC25A2, model="Blosum62")

# Write distance matrix to file
write.nexus.dist(dist_obj, "distance_matrix.csv", upper=FALSE, diag=FALSE)

# Make tree with Neighbor-Joining
tree <- NJ(dist_obj)

# Cladogram grouped by larger subtype
pdf("cladogram.pdf")
plt <- ggtree(tree, branch.length='none', layout='circular') +
       geom_tiplab(size=3) +
	   geom_highlight(node=62, fill="orange", alpha=.6) +
	   geom_highlight(node=79, fill="purple", alpha=.6) +
	   geom_highlight(node=78, fill="blue", alpha=.6) +
	   geom_highlight(node=72, fill="deeppink1", alpha=.6) +
	   geom_highlight(node=57, fill="deeppink2", alpha=.6) +
	   geom_highlight(node=65, fill="deeppink3", alpha=.6) +
#	   geom_text2(aes(subset=!istip, label=node), hjust=-.3, size=2) +
	   theme(plot.margin=unit(c(15,15,15,15), "mm"))
plt
dev.off()

# Cladogram grouped by more specific subtypes
pdf("SLC25A2_cladogram.pdf")
plt2 <- ggtree(tree, branch.length='none', layout='circular') +
        geom_tiplab(size=3) +
		geom_highlight(node=98, fill="darkgreen", alpha=.6) +
		geom_highlight(node=78, fill="green", alpha=.6) +
		geom_highlight(node=87, fill="blue", alpha=.6) +
		geom_highlight(node=96, fill="purple", alpha=.6) +
		geom_highlight(node=101, fill="hotpink", alpha=.6) +
		geom_highlight(node=90, fill="yellow", alpha=.6) +
		geom_highlight(node=102, fill="orange", alpha=.6) +
#		geom_text2(aes(subset=!istip, label=node), hjust=-.3, size=2) +
		theme(plot.margin=unit(c(15,15,15,15), "mm"))
plt2
dev.off()
