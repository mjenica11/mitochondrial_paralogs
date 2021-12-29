# Power calculations with PROPER
#BiocManager::install("PROPER")

library(PROPER)
# Simulate male and female counts from neg bin distribution using Gilad data
sim.opts.Gilad = RNAseq.SimOptions.2grp(ngenes=400, p.DE=0.05, lOD="gilad", 
										lBaselineExpr="cheung")

# Run simulations; evaluate power when there are Nreos samples in each 
# treatment group
simres = runSims(Nreps = c(3,4,6,10,11,13,18,35), sim.opts=sim.opts.Gilad, 
				 DEmethod="edgeR", nsims=100)

# Evaluate the powers
powers = comparePower(simres, alpha.type="fdr", alpha.nominal=0.1,
					  stratify.by="expr", delta=0.5)

summaryPower(powers)
plotPower(powers)
