library(ape)
library(phangorn)
library(seqinr)
library(ggtree)
library(ggrepel)

# Import data
carriers <- read.dna("~/mitochondrial_paralogs/named_SLC25_mRNA_mafft_alignment.clustalw", format="clustal")

# Convert to phyDat object
carriers_phyDat <- phyDat(carriers, type="DNA", levels=NULL)

# Replace SLC25 nomenclature with common names to make tree more human readable
common_names <- c("CIC", "DIC", "GC2", "GC1", "BAC", "HDCMP", "OGC", "APC2",
				  "APC3", "APC4", "CG169", "CoAPC", "SLC25A34", "SLC25A45",
				  "SLC25A35", "MFRN2", "SIDBA2", "SLC25A44", "UCP3", "UCP2",
				  "MTCH1", "ORC1", "ORC2", "TPC", "MFRN1", "AAC1", "AAC2",
				  "AAC3", "CAC", "SLC25A48", "AGC1", "AGC2", "UCP5", "APC1",
				  "GDA", "CFNC", "CG7943", "MCART2", "MCART6", "AAC4", "PNC2",
				  "MTCH2", "PCH1E", "UCP6", "MFT", "MCFP", "ODC", "UCP4", "PiC",
				  "PNC1", "SLC25A43", "SAMC", "UCP1")

names(carriers_phyDat) <- common_names

# Perform likelihood ratio test to decide on the model of nt evolution that
# best fits the data
#mt <- modelTest(carriers_phyDat)

# Highest logLik
#mt[which.max(mt$logLik),] # Model: GTR+G+I
#mt[which.min(mt$AIC),] # Model: GTR+G 
#mt[which.min(mt$BIC),] # Model: GTR+G 

# Make tree with general time reversible (GTR + G + I) model
#best_model <- mt$Model[which.min(mt$AICc)]
#env <- attr(mt, "env")
#fitStart <- eval(get("GTR+G+I", env), env)

# Plot tree created with GTR model
#plt <- plot(fitStart, use.edge.length=FALSE, main="mafft alignment with GTR+G+I model")

# Calculate distance matrix with JC69 model
dna_dist <- dist.ml(carriers_phyDat, model="JC69")

# Estimate trees using UPGMA and neighbor joining
carriers_UPGMA <- upgma(dna_dist)
carriers_NJ <- NJ(dna_dist)

# Which tree is a better fit for the data?
# Compare respective parsimony scores
parsimony(carriers_UPGMA, carriers_phyDat) # 55390
parsimony(carriers_NJ, carriers_phyDat) # 55123

# Replace SLC25 names with common names
carriers_UPGMA$tip.label <- common_names
carriers_NJ$tip.label <- common_names

# Plot trees
plot(carriers_UPGMA, main="alignment with mafft, JC69 substitution matrix, UPGMA clustering method")
plot(carriers_NJ, main="alignment with mafft, JC69 substitution matrix, neighbor joining clustering method")

# Make cladogram with ggtree package
tree <- ggtree(carriers_UPGMA, branch.length='none', layout='circular') 
tree + geom_text(aes(label=label), nudge_x=0.7, label.size=0.25)
