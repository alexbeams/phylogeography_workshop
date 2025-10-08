rm(list=ls())

library(diversitree)
library(ape)
library(phytools)

# The previous problem sets have illustrated the behavior of trees
# simulated from these models, some of their statistical behaviors,
# and now all that remains is to use fitted models to accomplish
# the goals of phylogeography, i.e. ancestral character estimation.

# Note: for most BiSSE users, the primary goal seems not to be ancestral
# reconstruction per se, but inferring relationships between states
# and diversification/extinction rates. As such, the documentation
# and support for the reconstruction methods that use BiSSE is more
# sparse than for ace or simmap.


# Let's start by simulating a tree and
# fit the model to it. We will then infer ancestral states and compare
# these with the known values from the simulation itself.

# simulate a BiSSE tree here

pars <- c(
        lambda0=1,
        lambda1=3,
        mu0=0.1,
        mu1=0.1,
        q01=0.1,
        q10=0.12)

tree <- tree.bisse(pars, max.taxa=100, x0=0)

par(mfrow=c(1,2))
plot(history.from.sim.discrete(tree, states=c(0,1)),
	tree, col=c('0'='black','1'='red'))


# the make.bisse function is a scalar function of two variables
#	of the form loglik = loglik(c(lambda,mu))

loglik <- make.bisse(tree, tree$tip.state)

loglikc <- constrain(loglik, lambda0~lambda1)

# What are the MLEs?
mles <- find.mle(loglik, pars)
mlesc <- find.mle(loglikc, pars[-2])


# these are the parmeter values that this infers:
mlepars <- mles$par

# We use the "marginal" reconstruction, which should be similar
# to the type of calculation used by ace:
asr_marginal <- asr.marginal(loglik, mlepars)

plot(tree, show.tip.label=F)
nodelabels(pie=t(asr_marginal), piecol=c("black", "red"), cex=0.7)
tip_matrix <- cbind(1-tree$tip.state, tree$tip.state)  # Convert to probabilities
tiplabels(pie=tip_matrix, piecol=c("black", "red"), cex=0.7)


# Exercises:
# 1. Simulate a collection of BiSSE trees (n=30) and use asr_marginal
#	to reconstruct ancestral states. For the same set of trees, use
#	ace as well, and compare the reconstructions to each other. Is it
#	obvious that marginal reconstruction with the BiSSE model out-
#	performs the simpler Markov models in ace? Are there different
#	combinations of parameters where that is the case? What if you
#	use trees with more taxa (say, n=100 or n=500?)

# 2. In your simulations, do you notice particular situations where
#	where reconstructions tend to be inaccurate (when does BiSSE
#	get "tricked")? Are there particular configurations of tip data
#	that present clear challenges for the BiSSE model?

# 3. Set diversification rates in both states equal to each other when simulating
#	yout trees, and likewise for the extinction rates. Compare models
#	that impose constraints to reduce the number of parameters to estimate
#	with models that impose no constraints, and explore the accuracy of
#	reconstructions in each case. How sensitive are reconstructions to
#	the parametric constraints you impose?




