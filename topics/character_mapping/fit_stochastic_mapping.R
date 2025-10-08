rm(list=ls())


library(coda) # for handling output from MCMC
library(phytools) # simmap functions, Mk models

# This code is complementary to the other one in
#	this subdirectory, and walks through how
#	to fit models using a character mapping
#	approach.

# Most of the phytools stuff is based on the Mk
#	model, which is a particular type of
#	Markov model. M stands for "Markov",
#	and "k" indicates the number of states.
#	(Lewis has an interesting discussion
#	relating to the use of conditional
#	likelihoods as a means of overcoming
#	acquisition bias -- the tendency to only
#	record changes of state, and not
#	maintenance of a character).


# fitMk fits the model to the data using nlminb,
#	and optimization routine similar to quasi-
#	Newton methods like BFGS. Compare this to ace.

# mcmcMk fits Mk models to trees using MCMC.


data(sunfish.tree)
data(sunfish.data)

# extract discrete character (fish diet):
sunfish_tipdata <-setNames(sunfish.data$feeding.mode,
    rownames(sunfish.data))

# fit "ER" model
fit.ER<-fitMk(sunfish.tree,sunfish_tipdata,model="ER")
print(fit.ER)


# This will run for a while, and you should see
#	traceplots updating in real time until it
#	finishes:
fit.mcmc <- mcmcMk(sunfish.tree, sunfish_tipdata,
	model="ER", ngen=10000)


# the dimensions of fit.mcmc will be
#	ngen * (number of parameters in the Mk model + 2),
# 	which for this run will be 10000 x 3 because
#	we are only estimating a single parameter

# The first column of fit.mcmc is timesteps of the
#	MCMC routine. The last column of fit.mcmc is
#	the logLikelihood, and the columns between first
#	and last are the parameters at each timestep.

# We can look at the traceplots:
plot(fit.mcmc)

# We can also plot the posterior distribution
#	for our parameters. By default, this should
#	toss out the first 20% of the MCMC updates
# 	as part of a burn-in (we discard transients
#	because we only care about the stationary
#	distribution of the MCMC run)

plot(density(fit.mcmc))

# Of note, we can also use base R plotting to accomplish
#	these things as well. Here is a histogram of
#	the parameter values of the entire MCMC run:

hist(fit.mcmc[,'[pisc,non]'],breaks=100)

# and the last 80% of the run:
hist(fit.mcmc[c(2001:10000),'[pisc,non]'],breaks=100)

# I also really like to plot profile log-
#	likelihood (logLik~parameters):

plot(fit.mcmc[,'logLik']~fit.mcmc[,'[pisc,non]'])

# ^^ This helps us visualize the peak of the
#	likelihood surface. Can be useful for
#	diagnosing convergence problems

####################################################
####################################################
####################################################

# Exercise set 1:
#	1. Did the fitMk routine produce the same
#	estimates as mcmcMk? How would you determine
#	this?

#	2. Running mcmcMk gives you distributions
#	over the parameters. This is nice for
#	visualizing/quantifying uncertainty in the
#	estimates. Does fitMk produce something
#	analagous?

#	3. What was assumed about the root states in
#	these models? Examine the help pages for fitMk
#	and mcmcMk, and re-run results for a different
#	choice of the "root prior". Are estimates of
#	the transition rate sensitive to your choice?

####################################################
####################################################
####################################################


# We can fit a more complicated model with two rates:

fit.mcmc_2rates <- mcmcMk(sunfish.tree, sunfish_tipdata,
	model="ARD", ngen=10000)

# Now, we have profile log-likelihoods in two parameters:

par(mfrow=c(2,1))
plot(fit.mcmc_2rates[,'logLik']~
	fit.mcmc_2rates[,'[pisc,non]'])
plot(fit.mcmc_2rates[,'logLik']~
	fit.mcmc_2rates[,'[non,pisc]'])


# We can also visualize correlations in parameters:
par(mfrow=c(1,1))
plot(fit.mcmc_2rates[,'[pisc,non]'],
	fit.mcmc_2rates[,'[non,pisc]'])



# We're going to write down a function that accepts the two
#	transition rates and spits out a character mapping
#	(a make.simmap output).
#	We'll want to estimate the root state probabilities too,
#	but for now we can just set them to something fixed.

# This is just so to make the algorithm below a little easier
#	to read.

get_simmap <- function(q12, q21, tree, tipdata){

	Q <- matrix(c(-q12,q12,
		      q21,-q21),2,2)
	rownames(Q) <- colnames(Q) <- levels(tipdata)

	# We simulate from the root to the tips.
	#	We have to specify an initial condition
	#	at the root:

	Pi <- c(0,1)

	# let's just simulate 1 character mapping for
	#	each parameter update for now; can
	#	change this later

	sims <- make.simmap(tree, tipdata, Q=Q, pi=Pi, nsim=1)
	return(sims)
}

# let's now simulate a mapping for each step of our MCMC:
# (but let's just use the last 1000 rows of the fit.mcmc_2rates
#	object)

simmaps <- apply(fit.mcmc_2rates[9001:10000,], 1, function(x) {
	get_simmap(x[2],x[3], sunfish.tree, sunfish_tipdata)

})

# Let's start by summarizing the number of mutations that happen across the whole tree:
counts <- lapply(simmaps, countSimmap)
numMutations <- sapply(counts, function(x) x$N)
transitions <- lapply(counts, function(x)x$Tr)
transitions_non_pisc <- sapply(transitions, function(x) x[1,2])
transitions_pisc_non <- sapply(transitions, function(x) x[2,1])

### plot an example character mapped tree, and histograms of number of mutations
sim = simmaps[[1]]
#  overall, and by type
par(mfrow=c(2,2))
## plot the character mapped tree
sim$tip.label <- paste(sim$tip.label, "(", sunfish_tipdata, ")", sep="")
plotSimmap(sim, fsize = 0.8)
## plottig make.simmap objects messes up margins, let's reset these:
par(mar=c(5,4,4,1))
## histogram of total mutations:
hist(numMutations)
## histogram of mutations non -> pisc
hist(transitions_non_pisc)
## histogram of mutations pisc -> non
hist(transitions_pisc_non)

####################################################
####################################################
####################################################
# Exercise set 2
####################################################
####################################################
####################################################

# 1. Why do you think the mcmcMk routine does not
#	propose character mappings alongside parameters?
#	Would anything be gained by treating character
#	mappings as latent variables in MCMC
#	and updating those alongside parameters?

# Check this paper out if you are interested in MCMC
#	and character mapping:
#
#  Mutations as Missing Data: Inferences on the Ages and
#	Distributions of Nonsynonymous and Synonymous
#	Mutations,
#	by Rasmus Nielsen
#	Genetics 159: 401-411 (September 2001)


# 2. Comment on any interesting things you see in profile
#	loglikelihoods for models with multiple parameters.
#	Can you deduce anything interesting about the geometry
#	of the likelihood surface?

# 3. In the previous runs, how did we estimate Pi, the
#	probability of the root states? Examine the
#	fitMk and mcmcMk help page to determine what the
#	default is. Try fitting the model under a different
#	setting to see how results change.

### Add in trees with more tipstates as exercises.
