##==============================================
## Sample Size Estimation
##==============================================
##
## Starting off with r = 10, 100, 1000
## H0: a = 50
##     b = 50
## H1: a = 49, 47, 45
##     b = 1, 3, 5
##==============================================
## Author: samarov
###############################################################################
require("rjags")
require("multicore")

sampBin <- function(r = NULL, epsilon = NULL,
		a = NULL, b = NULL){
	E <- rbinom(n=r, size=1, prob=epsilon)
	B <- rep(x=NA, times=r)
	p <- a/(a+b)
	for (i in 1:r) {
		
		if (E[i] == 0) {
			B[i] <- rbinom(n=1, size=1, prob=p)
		} else {
			B[i] <- rbinom(n=1, size=1, prob=(1 - p))
		}
	}
	return(list(B = B, E = E))
}


sampleDist <- function(r, a, b, a1 = NULL, b1 = NULL,
		epsilon, iter = 100, flat = FALSE, spike = 1){
	
	model <- "
			model{
			for(i in 1:r) {
			prob[i] <- (1 - E[i])*p.b + E[i]*(1 - p.b)
			B[i] ~ dbern(prob[i])
			E[i] ~ dbern(epsilon)
			}
			p.b ~ dbeta(a, b)
			}
			"
	
	tmpf <- tempfile()
	tmps <- file(tmpf, "w")
	cat(model, file = tmps)
	close(tmps)
	
	if(is.null(a1) | is.null(b1)){
		a1 <- a; b1 <- b
	}
	
	ar <- a * spike
	br <- spike - ar
	a1r <- a1 * spike
	b1r <- spike - a1r
	
	iterIndx <- cbind(1:iter)
	
	st <- system.time({simRes <- multicore::mclapply(iterIndx, mc.cores = 16, FUN = function(u){
							
							sb1 <- sampBin(r = r, epsilon = epsilon,
									a = a, b = b)
							sb2 <- sampBin(r = r, epsilon = epsilon,
									a = a1, b = b1)
							
							if(!flat){
								as1 <- sum(sb1$B)
								bs1 <- ifelse(as1 == r, 1, r - as1)
								
								as2 <- sum(sb2$B)
								bs2 <- ifelse(as2 == r, 1, r - as2)
								
								as1 <- spike * as1; bs1 <- spike * bs1
								as2 <- spike * as2; bs2 <- spike * bs2
								
								inits1 <- function() list(p.b = rbeta(1, as1, bs1))
								inits2 <- function() list(p.b = rbeta(1, bs2, bs2))
								
								data1 <- list(B = sb1$B, r = r, epsilon = epsilon, a = as1,
										b = bs1)
								data2 <- list(B = sb2$B, r = r, epsilon = epsilon, a = as2,
										b = bs2)
							}
							else{
								inits1 <- function() list(p.b = rbeta(1, 1, 1))
								inits2 <- function() list(p.b = rbeta(1, 1, 1))
								data1 <- list(B = sb1$B, r = r, epsilon = epsilon, a = 1,
										b = 1)
								data2 <- list(B = sb2$B, r = r, epsilon = epsilon, a = 1,
										b = 1)
							}
							
							mod1 <- jags.model(tmpf, data = data1, inits = inits1, n.chains = 3, quiet = TRUE, n.adapt = 100)
							update(mod1, 1000)
							pSim1 <- coda.samples(mod1, "p.b", 1000)
							
							mod2 <- jags.model(tmpf, data = data2, inits = inits2, n.chains = 3, quiet = TRUE, n.adapt = 100)
							update(mod2, 1000)
							pSim2 <- coda.samples(mod2, "p.b", 1000)
							
							list(sampEst1 = pSim1, sampEst2 = pSim2)
							
						})})
	
	attr(simRes, "time") <- st
	
	return(simRes)
	
}


powerSim <- function(a0, b0, a1, b1, r, epsilon, iter, flat, spike){
	
	nn <- length(a1)
	nr <- length(r)
	ne <- length(epsilon)
	
	res <- vector("list", ne)
	
	for(i in 1:ne){
		res[[i]] <- vector("list", nr)
		names(res[[i]]) <- paste("read depth", r)
		for(j in 1:nr){
			res[[i]][[j]] <- vector("list", nn + 1)
			names(res[[i]][[j]]) <- c("H0", paste("H1", 1:nn))
		}
	}
	names(res) <- paste("epsilon", epsilon)
	
	for(k in 1:ne){
		for(i in 1:nr){
			s0 <- sampleDist(r = r[i], a = a0, b = b0,
					epsilon = epsilon[k], iter = iter, flat = flat,
					spike = spike)
			res[[k]][[i]][[1]] <- s0
			for(j in 1:nn){
				s1 <- sampleDist(r = r[i], a = a0, b = b0,
						a1 = a1[j], b1 = b1[j],
						epsilon = epsilon[k], iter = iter, flat = flat,
						spike = spike)
				res[[k]][[i]][[(j + 1)]] <- s1
			}
		}
	}
	
	return(res)
	
}

## Model parameters
a0 <- .5
b0 <- .5
a1 <- c(48, 46, 44, 42, 40, 30, 20, 10)/100
b1 <- 1 - a1
r <- c(10, 30, 100, 200)
epsilon <- c(.001, .01)

## Simulating data for power calculation estimates
## Weak
stWeak <- system.time({powerSim5050Weak <- powerSim(a0 = a0, b0 = b0, a1 = a1, b1 = b1,
					r = r, epsilon = epsilon, iter = 500, flat = TRUE, spike = 1)})
save(powerSim5050Weak, file = "output/powerSim5050Weak.RData")

## Medium
stMed <- system.time({powerSim5050Med <- powerSim(a0 = a0, b0 = b0, a1 = a1, b1 = b1,
					r = r, epsilon = epsilon, iter = 500, flat = FALSE, spike = 1)})
save(powerSim5050Med, file = "output/powerSim5050Med.RData")
