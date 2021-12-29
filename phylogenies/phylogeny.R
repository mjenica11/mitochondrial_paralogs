library(ape)
library(phangorn)
library(seqinr)
library(ggtree)

# Import data
carriers <- read.dna("~/Downloads/alignments/all_carriers.fasta", format="fasta")

# Convert to phyDat object
carriers_phyDat <- phyDat(carriers, type="DNA", levels=NULL)

# Perform likelihood ratio test to decide on the model of nt evolution that
# best fits the data
mt <- modelTest(carriers_phyDat)

# Highest logLik
mt[which.max(mt$logLik),] # Model: GTR+I
mt[which.min(mt$AIC),] # Model: HKY+G
mt[which.min(mt$BIC),] # Model: HKY+G

# Couldn't figure out which function models distance with GTR so using
# Jukes Cantor instead for now...

# Calculate distance matrix with JC69 model
dna_dist <- dist.ml(carriers_phyDat, model="JC69")

# Estaimate trees using UPGMA and neighbor joining
carriers_UPGMA <- upgma(dna_dist)
carriers_NJ <- NJ(dna_dist)

# Plot trees
plot(carriers_UPGMA, main="UPGMA")
plot(carriers_NJ, main="Neighbor Joining")

# Which tree is a better fit for the data?
# Compare respective parsimony scores
parsimony(carriers_UPGMA, carriers_phyDat) # 69546
parsimony(carriers_NJ, carriers_phyDat) # 69543

# Make cladogram with ggtree package
tree <- ggtree(carriers_UPGMA, branch.length='none', layout='circular') + geom_tiplab()
tree


