rm(list=ls())

library(diversitree)

# State-dependent speciation and extion (SSE) models
#	are the subject of this set of exercises. The
#	workhorse of these methods are the various flavors
#	of birth-death models for generating trees.

# So far, we have been treating the trees as given 
#	without asking how they are created (by some
#	biologist, presumably). If the traits we are
#	interested in do not influence diversification, 
#	extinction/death, or the probability we observe
#	a particular lineage, then the standard Markov models
#	on trees seem sufficient. In some cases, traits
#	might influence the process creating the phylogenetic
#	relationships we observe, and we would like to detect
#	these effects.

# It is for this reason that we start by considering models for
#	tree generation (sometimes called "tree priors" by
#	people who use BEAST, which we will get to later in 
#	the week).

# Let's start by simulating some trees using birth-death models.

tree <- tree.bd(c(lambda=10, mu=5), max.taxa = 50)
plot(tree)

# Neat-o. But it doesn't always work. Let's use the lapply function
#	to generate a list of simulated outputs (like a for loop,
#	but better).

numtrees <- 50
trees <- lapply(1:numtrees, 
		function(x) tree.bd(c(lambda=10, mu=5), max.taxa = 50)
)

dieoffs <- sapply(trees, is.null)
numdieoffs <- sum(dieoffs)

# Sometimes we fail to generate a tree with max.taxa number of 
#	lineages because they all go extinct (this is stochastic).
#	Let's see how many we were able to generate:

# Exercises:
# 1. The parameter "lambda" is the birth/diversification rate,
#	and the parameter "mu" is the mortality/extinction rate.
#	Try different combinations of the parameters (try 
#	lambda = 20, 50, 100, and mu = 20, 50, 100). What do you
#	notice about your simulations?
#	Hint: you should write a function called "getnumdieoffs"
#	that accepts lambda, mu, and numtaxa, and then spits out 
#	trees. Then you can create another function that acts on 
#	trees to return numdieoffs.
#	Organize the results of your simulations in a nice display.




