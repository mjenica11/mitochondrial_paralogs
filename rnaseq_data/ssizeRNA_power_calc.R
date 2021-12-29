# Calculate sample size needed

# Install ssizeRNA package and load library
#BiocManager::install('ssizeRNA')
library(ssizeRNA)

SAMPLE_SIZE <- 19
# Only parameter modified: m = sample size per treatment group
# Benjamini and Hochberg method
set.seed(711)
check.power(nGenes = 200, m = SAMPLE_SIZE, mu = 10, disp = 0.1, fc = 2, sims = 100)

# Sample size calculations for two-sample RNA-seq experiments with a 
# single set of parameters
set.seed(711)
size1 <- ssizeRNA_single(nGenes = 70, pi0 = 0.8, m = 100, mu = 10, disp = 0.1,
						 fc = 2, fdr = 0.5, power = 0.8, maxN = SAMPLE_SIZE)
size1$power
