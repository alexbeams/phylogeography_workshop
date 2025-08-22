rm(list=ls())

library(phytools)

#set.seed(123)

#no. of stochastic character mappings to simulate:
nsim = 2000

tree <- pbtree(n=40, scale=1)

# Transition rate matrix (two states A,B)
Q <- matrix(c(-1,1,
              1,-1),2,2)
rownames(Q) <- colnames(Q) <- c("A","B")

# can scale down Q to have fewer transitions:
Q <- Q/1

# Function to simulate character data and run SCM for one tree
tip_states <- sim.Mk(tree, Q)
sims <- make.simmap(tree, tip_states, model="ER", nsim=nsim)

#for(i in 1:nsim){
#	plot(sims[[i]])
	Sys.sleep(.002)
#}

# Let's start by summarizing the number of mutations that happen across the whole tree:

counts <- lapply(sims, countSimmap)

numMutations <- sapply(counts, function(x) x$N)

transitions <- lapply(counts, function(x)x$Tr)

transitions_AB <- sapply(transitions, function(x) x[1,2])

transitions_BA <- sapply(transitions, function(x) x[2,1])

par(mfrow=c(2,2))

hist(numMutations)

hist(transitions_AB)

hist(transitions_BA)

plot(sims[[1]])

