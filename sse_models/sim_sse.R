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

tree <- tree.bd(c(lambda=1, mu=0), max.taxa = 50)
plot(tree)

# Neat-o. But it doesn't always work if mu>0. 
#	Let's use the lapply function to generate a 
#	list of simulated outputs (lapply is like a 
#	for loop, but more high-brow).

numtrees <- 50
treelist <- lapply(1:numtrees, 
		function(x) tree.bd(c(lambda=10, mu=5), max.taxa = 50)
)

dieoffs <- sapply(treelist, is.null)
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
#	treelist. Then you can create another function that acts on 
#	treelist to return numdieoffs.
#	Organize the results of your simulations in a nice display.

####
# Solution
###

lambdas <- c(20, 50, 100)
mus <- c(20, 50, 100)
numtaxas <- c(10, 20, 30, 40, 50)


gettrees <- function(lambda, mu, numtaxa, numtrees){
	treelist <- lapply(1:numtrees, 
		function(x){ 
			tree.bd(c(lambda=lambda, mu=mu), 
				max.taxa = numtaxa)
		})
	return(treelist)
}

pars <- expand.grid(lambdas, mus, numtaxas)
colnames(pars) <- c('lambda','mu','numtaxa')

getnumdieoffs <- function(lambda, mu, numtaxa, numtrees){
	treelist <- gettrees(lambda, mu, numtaxa, numtrees)
	numdieoffs <- sum(sapply(treelist, is.null))
	return(numdieoffs)
}

# let's count number of dieoffs for 100-taxa trees:
simdieoffs <- cbind(pars, apply(pars, 1, function(x){
	getnumdieoffs(x[1],x[2],x[3],100)
}))

colnames(simdieoffs)[3] <- c('numtaxa')
colnames(simdieoffs)[4] <- c('numdieoffs')


# 2. Try to simulate trees with lambda = 1, mu = 2, and 10 taxa.
#	How many simulations do you need to run until you start
#	seeing trees with 10 taxa produced? 


# 3. Use ?trees to examine the help page for some of other tree
#	simulation models. You will notcie tree.bd and tree.yule
#	toward the bottom. Use your investigative abilities to 
#	research the Yule process and compare it to the birth-
#	death process. 

# 4. Do tree.yule(1, max.taxa=10) and tree.bd(lambda=2,mu=1, max.taxa=10)
#	simulate the same process? Justify your answer.

# 5. Do tree.yule(1, max.taxa=10) and tree.bd(lambda=1,mu=0, max.taxa=10)
#	simualte the same process? Justify your answer.

# 6. The diversitree package can simulate tip data using Markov
#	models for two discrete states initialized from a particular 
#	state at the root. One model is the "mk2" model: 
#	sim.character(tree, pars, x0=0, model="mk2", br=NULL)
#	the "pars" argument is a vector of the form (q12,q21)
#	corresponding to the transition rates of the stochastic
#	rate matrix, Q (the model argument can also be set to "mkn"
#	to simulate Markov models with n states).
#	Compare sim.character with the simMk function from phytools. 
#	Confirm whether these two R functions simulate the same process.	

## Add some exercises using max.t instead of max.taxa


#########
# BISSE
#########

# So now, we are able to simulate trees using the birth-death process,
#	and we learned earlier this week how to simulate trait evolution
#	on a given tree. Now, we want to explore what happens if 
#	the tree generation process interacts with the traits. To do this,
#	we will simulate Binary State Speciation and Extinction (BiSSE)
#	models.

# BiSSE has 5 parameters: birth/diversification rates (lambda0, lambda1)
#	corresponding to the two values of the trait (0 and 1); two 
#	death/extinction rates (mu0, mu1), and rates of the stochastic
#	rate matrix, Q, (q01 and q10).

pars <- c(
	lambda0=0.1, 
	lambda1=0.2, 
	mu0=0.03, 
	mu1=0.03, 
	q01=0.01, 
	q10=0.01) 

tree <- tree.bisse(pars, max.taxa=50, x0=0)

# State change information is available in tree$node.state,
#	tree$tip.state, and tree$edge.state, but for plotting
#	it is convenient to invoke the following:

plot(history.from.sim.discrete(tree, states=c(0,1)),
	tree, col=c('0'='black','1'='red')) 

# I'm going to write a function to simulate and plot 
#	from parameter values I supply:


## to do: need to figure out why the colors don't seem to match the states:

plotbissesim <- function(
	lambda0=1,
	lambda1=1,
	mu0=0,
	mu1=0,
	q01=1,
	q10=1,
	maxtaxa=50,
	x0=0){

	pars <- c(
		lambda0=lambda0, 
		lambda1=lambda1, 
		mu0=mu0, 
		mu1=mu1, 
		q01=q01, 
		q10=q10) 

	tree <- tree.bisse(pars, max.taxa=maxtaxa, x0=x0)
	plot(history.from.sim.discrete(tree, states=c(0,1)),
		tree, col=c('0'='black','1'='red')) 
	legend("topleft", 
           legend=c("state 0","state 1"), 
           col=c("black","red"), 
           lty=1, lwd=2, 
           bty="n", cex=0.8)
}

# Because I supplied defualt values, this will just run:
plotbissesim()

# Or, I can tweak parameters (as many at a time as I want):
plotbissesim(lambda0=10)

plotbissesim(lambda0=0.1,lambda1=1,q01=.1,q10=1,x0=1)

# BiSSE Simulation Exercises:

# Simulate a population that starts in a hostile environment but can move
#	to a safer habitat with less predation. Assume the organisms are dumb 
#	(or just plants), so that individuals tend to spend the same amount of
#	time in each habitat before moving to the other one. Do the simulations 
#	match your intuition about what should happen to the population in 
#	this scenario? Do you see different outcomes if migration is slow
#	or fast relative to the birth and death rates? What happens if you 
#	simulate more taxa? (Try up to 10^4 or 10^5).

# Simulate a population that can move to a new habitat where there are 
#	more resources available, and organisms can produce more offspring.
#	Assume that organisms have the same life expectancy in each area.
#	You can also assume that individuals spend the same amount of
#	time in each location. How does this scenario compare to the previous
#	one, in terms of the trees you generate? Try different numbers of taxa.

# Now, suppose through no fault of their own that organisms are transported
#	to a terrible environment at a much faster rate than they can
#	escape. If you didn't know the parameters you used, would it be
#	obvious from the trees that one of the locations is worse than the
#	the other?
 
# Supose the organisms are intelligent, and preferentially move to the 
#	location where there is less predation, and/or better resources.
#	Do the phylogenies contain information about the crummy habitat?



