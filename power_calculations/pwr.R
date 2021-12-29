library(pwr)

# Power calculations for balanced one-way ANOVA tests
# Purpose: Compute power of test or determine parameters to obtain target 
# power
# Parameters:
# k = number of groups
# n = number of observations per group
# f = effect size
# sig.level = Type I error probability
# power = 1 minus Type II error of probability

SAMPLE_SIZE <- 13
N_GENES <- 200

pwr.anova.test(k = N_GENES, n = SAMPLE_SIZE, f = 2, sig.level = 0.05,
			   power = NULL)

