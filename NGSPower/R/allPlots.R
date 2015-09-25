# Aggregating plots
# 
# Author: samarov
###############################################################################
library(ggplot2)
library(ks)

##=======================================================================
## Analysis of 99.99% purity with weak prior
##=======================================================================
source("R/analysis_99_1_weak.R")
## Densities
png("plots/eps001_99pure_weak.png", height = 800, width = 800)
pEps001 <- ggplot(sdiffsDf, aes(x = Difference, y = Density, colour = Hypothesis))
pEps001 + geom_line() + facet_grid(Reads~., scales = "free", space = "free_y") + 
		opts(title = "Flat prior distribution")
dev.off()
## Power plots
png("plots/pow_99pure_weak.png", height = 800, width = 800)
powWeak <- ggplot(powDf, aes(x = `Sample 2 Purity`, y = Power, colour = `Read Depth`))
powWeak + geom_line() + facet_grid(`ErrorRate`~.) + 
		opts(title = "Flat prior distribution")
dev.off()

##=======================================================================
## Analysis of 99.99% purity with moderate prior
##=======================================================================
source("R/analysis_99_1_mod.R")
## Densities
png("plots/eps001_99pure_mod.png", height = 800, width = 800)
pEps001 <- ggplot(sdiffsDf, aes(x = Difference, y = Density, colour = Hypothesis))
pEps001 + geom_line() + facet_grid(Reads~., scales = "free", space = "free_y") + 
		opts(title = "Moderate prior distribution")
dev.off()
## Power plots
png("plots/pow_99pure_mod.png", height = 800, width = 800)
powMod <- ggplot(powDf, aes(x = `Sample 2 Purity`, y = Power, colour = `Read Depth`))
powMod + geom_line() + facet_grid(`ErrorRate`~.) + 
		opts(title = "Moderate prior distribution")
dev.off()

##=======================================================================
## Analysis of 50% purity with weak prior
##=======================================================================
source("R/analysis_50_50_weak.R")
## Densities
png("plots/eps001_50pure_weak.png", height = 800, width = 800)
pEps001 <- ggplot(sdiffsDf, aes(x = Difference, y = Density, colour = Hypothesis))
pEps001 + geom_line() + facet_grid(Reads~., scales = "free", space = "free_y") + 
		opts(title = "Flat prior distribution")
dev.off()
## Power plots
png("plots/pow_50pure_weak.png", height = 800, width = 800)
powMod <- ggplot(powDf, aes(x = `Sample 2 Purity`, y = Power, colour = `Read Depth`))
powMod + geom_line() + facet_grid(`ErrorRate`~.) + 
		opts(title = "Moderate prior distribution")
dev.off()

##=======================================================================
## Analysis of 50% purity with moderate prior
##=======================================================================
source("R/analysis_50_50_mod.R")
## Densities
png("plots/eps001_50pure_mod.png", height = 800, width = 800)
pEps001 <- ggplot(sdiffsDf, aes(x = Difference, y = Density, colour = Hypothesis))
pEps001 + geom_line() + facet_grid(Reads~., scales = "free", space = "free_y") + 
		opts(title = "Moderate prior distribution")
dev.off()
## Power plots
png("plots/pow_50pure_mod.png", height = 800, width = 800)
powMod <- ggplot(powDf, aes(x = `Sample 2 Purity`, y = Power, colour = `Read Depth`))
powMod + geom_line() + facet_grid(`ErrorRate`~.) + 
		opts(title = "Moderate prior distribution")
dev.off()