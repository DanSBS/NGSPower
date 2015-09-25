# Analysis file for the 99 - 1 data
# 
# Author: samarov
###############################################################################
library(ggplot2)
library(ks)
## Load simulation results for weak assumptions
load("output/powerSim5050Weak.RData")

## Model parameters
a0 <- .5
b0 <- .5
a1 <- c(48, 46, 44, 42, 40, 30, 20, 10)/100
b1 <- 1 - a1
r <- c(10, 30, 100, 200)
epsilon <- c(.001, .01)
ne <- length(epsilon)
nt <- length(a1) + 1

## Assuming there is approximately a .05 difference between samples
fdrAdj <- 1e-7
## Level of significance (looking at one sided tests here)
alpSig <- 0.05
## Getting the differences of the null and alternative hypotheses
diffs10 <- lapply(powerSim5050Weak, function(u){
			lapply(u[[1]], function(v){
						
						unlist(lapply(v, function(w){
											if(length(w) > 1) 
												unlist(w[[1]]) - unlist(w[[2]])
										}))
					})
		})


pow10 <- lapply(diffs10, function(u){
			qalpha <- quantile(u[[1]], 1 - alpSig * fdrAdj)
			out <- rep(NA, length(u) - 1)
			for(j in 2:length(u)){
				out[j-1] <- sum(u[[j]] >= qalpha)/length(u[[j]])
			}
			out
		})

diffs30 <- lapply(powerSim5050Weak, function(u){
			lapply(u[[2]], function(v){
						unlist(lapply(v, function(w){
											if(length(w) > 1) 
												unlist(w[[1]]) - unlist(w[[2]])
										}))
					})
		})


pow30 <- lapply(diffs30, function(u){
			qalpha <- quantile(u[[1]], 1 - alpSig * fdrAdj)
			out <- rep(NA, length(u) - 1)
			for(j in 2:length(u)){
				out[j-1] <- sum(u[[j]] >= qalpha)/length(u[[j]])
			}
			out
		})


diffs100 <- lapply(powerSim5050Weak, function(u){
			lapply(u[[3]], function(v){
						unlist(lapply(v, function(w){
											if(length(w) > 1) 
												unlist(w[[1]]) - unlist(w[[2]])
										}))
					})
		})

pow100 <- lapply(diffs100, function(u){
			qalpha <- quantile(u[[1]], 1 - alpSig * fdrAdj)
			out <- rep(NA, length(u) - 1)
			for(j in 2:length(u)){
				out[j-1] <- sum(u[[j]] >= qalpha)/length(u[[j]])
			}
			out
		})

diffs200 <- lapply(powerSim5050Weak, function(u){
			lapply(u[[4]], function(v){
						unlist(lapply(v, function(w){
											if(length(w) > 1) 
												unlist(w[[1]]) - unlist(w[[2]])
										}))
					})
		})
pow200 <- lapply(diffs200, function(u){
			qalpha <- quantile(u[[1]], 1 - alpSig * fdrAdj)
			out <- rep(NA, length(u) - 1)
			for(j in 2:length(u)){
				out[j-1] <- sum(u[[j]] >= qalpha)/length(u[[j]])
			}
			out
		})
## Looking at the smoothed histogram
##
sdiffs10 <- lapply(diffs10, function(u){
			lapply(u, function(v) {
						if(!is.null(v)) 
							bkde(v, range.x = range(v))
						else
							NULL
					})
		})
xrng10 <- lapply(sdiffs10, function(u){
			range(unlist(lapply(u, function(v) v$x)))
		})
yrng10 <- lapply(sdiffs10, function(u){
			range(unlist(lapply(u, function(v) v$y)))
		})

##
sdiffs30 <- lapply(diffs30, function(u){
			lapply(u, function(v) bkde(v, range.x = range(v)))
		})
xrng30 <- lapply(sdiffs30, function(u){
			range(unlist(lapply(u, function(v) v$x)))
		})
yrng30 <- lapply(sdiffs30, function(u){
			range(unlist(lapply(u, function(v) v$y)))
		})

sdiffs100 <- lapply(diffs100, function(u){
			lapply(u, function(v) bkde(v, range.x = range(v)))
		})
xrng100 <- lapply(sdiffs100, function(u){
			range(unlist(lapply(u, function(v) v$x)))
		})
yrng100 <- lapply(sdiffs100, function(u){
			range(unlist(lapply(u, function(v) v$y)))
		})

sdiffs200 <- lapply(diffs200, function(u){
			lapply(u, function(v) bkde(v, range.x = range(v)))
		})
xrng200 <- lapply(sdiffs200, function(u){
			range(unlist(lapply(u, function(v) v$x)))
		})
yrng200 <- lapply(sdiffs200, function(u){
			range(unlist(lapply(u, function(v) v$y)))
		})

## Placing these into data frames for visualizations using ggplot2
sdiffs10Df <- rbind()
for(i in 1:length(sdiffs10[[1]])){
	if(i == 1)
		hyp <- "99.99%"
	else
		hyp <- paste(100* a1[i-1], "%", sep = "")
	
	sdiffs10Df <- rbind(sdiffs10Df, 
			data.frame(Difference=sdiffs10[[1]][[i]]$x, 
					Density=sdiffs10[[1]][[i]]$y, Hypothesis=hyp,Epsilon="0.001",
					Reads="10"))
}

sdiffs30Df <- rbind()
for(i in 1:length(sdiffs30[[1]])){
	if(i == 1)
		hyp <- "99.99%"
	else
		hyp <- paste(100* a1[i-1], "%", sep = "")
	
	sdiffs30Df <- rbind(sdiffs30Df, 
			data.frame(Difference=sdiffs30[[1]][[i]]$x, 
					Density=sdiffs30[[1]][[i]]$y, Hypothesis=hyp,Epsilon="0.001",
					Reads="30"))
}

sdiffs100Df <- rbind()
for(i in 1:length(sdiffs100[[1]])){
	if(i == 1)
		hyp <- "99.99%"
	else
		hyp <- paste(100* a1[i-1], "%", sep = "")
	
	sdiffs100Df <- rbind(sdiffs100Df, 
			data.frame(Difference=sdiffs100[[1]][[i]]$x, 
					Density=sdiffs100[[1]][[i]]$y, Hypothesis=hyp,Epsilon="0.001",
					Reads="100"))
}


sdiffs200Df <- rbind()
for(i in 1:length(sdiffs200[[1]])){
	if(i == 1)
		hyp <- "99.99%"
	else
		hyp <- paste(100* a1[i-1], "%", sep = "")
	
	sdiffs200Df <- rbind(sdiffs200Df, 
			data.frame(Difference=sdiffs200[[1]][[i]]$x, 
					Density=sdiffs200[[1]][[i]]$y, Hypothesis=hyp,Epsilon="0.001",
					Reads="200"))
}

sdiffsDf <- rbind(sdiffs10Df, sdiffs30Df, sdiffs100Df, sdiffs200Df)
names(sdiffsDf)
# [1] "Difference" "Density"    "Hypothesis" "Epsilon"    "Reads"     
sdiffsDf$Reads <- paste("r =", sdiffsDf$Reads)

sdiffsDf$Hypothesis <- as.factor(sdiffsDf$Hypothesis)
sdiffsDf$Reads <- factor(sdiffsDf$Reads, levels = c("r = 10", "r = 30", "r = 100", "r = 200"))
sdiffsDf$Epsilon <- as.factor(sdiffsDf$Epsilon)

## Looking at plots of the power calcs
powEps001 <- c(pow10[[1]], pow30[[1]], pow100[[1]], pow200[[1]])
powEps01 <- c(pow10[[2]], pow30[[2]], pow100[[2]], pow200[[2]])
sampS <- paste("r =", rep(c(10, 30, 100, 200), each = length(pow10[[1]])))


powDf <- data.frame(Power = c(powEps001, powEps01), ss = sampS, 
		er = rep(c("0.001", "0.01"), each = length(powEps001)),
		ap = a1)
names(powDf)[2:4] <- c("Read Depth", "ErrorRate", "Sample 2 Purity")
powDf$`Read Depth` <- factor(powDf$`Read Depth`, levels = c("r = 10", "r = 30", "r = 100", "r = 200"))
powDf$`ErrorRate` <- as.factor(powDf$`ErrorRate`)
